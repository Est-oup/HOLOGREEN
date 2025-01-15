import os
import sys
import argparse
from Bio import SeqIO
import pandas as pd
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Assembler des contigs par bin avec Flye et organiser les sorties."
    )
    parser.add_argument("--fasta", required=True, help="Fichier FASTA contenant les contigs.")
    parser.add_argument("--clusters", required=True, help="Fichier TSV contenant les clusters et contigs (cluster.tsv).")
    parser.add_argument("--output", required=True, help="Dossier de sortie principal.")
    parser.add_argument("--flye_threads", default=1, type=int, help="Nombre de threads pour Flye (défaut: 1).")
    parser.add_argument("--max_parallel", default=1, type=int, help="Nombre maximum de bins à assembler en parallèle (défaut: 1).")
    return parser.parse_args()

def setup_output_dirs(output_dir):
    bins_dir = os.path.join(output_dir, "bins")
    flye_dir = os.path.join(output_dir, "flye_results")
    final_dir = os.path.join(output_dir, "final_output")
    logs_dir = os.path.join(output_dir, "logs")
    
    for dir_path in [bins_dir, flye_dir, final_dir, logs_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    return bins_dir, flye_dir, final_dir, logs_dir

def setup_logging(log_dir):
    log_file = os.path.join(log_dir, "process.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return log_file

def extract_contigs_by_bin(fasta_file, cluster_file, bins_dir):
    logging.info("Début de l'extraction des contigs par bin.")
    clusters = pd.read_csv(cluster_file, sep="\t", header=None, names=["bin", "contig"])
    fasta_sequences = {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}
    bins = clusters.groupby("bin")

    bin_fasta_files = []
    single_contig_bins = []
    for bin_id, group in bins:
        bin_fasta_file = os.path.join(bins_dir, f"bin_{bin_id}.fasta")
        with open(bin_fasta_file, "w") as f_out:
            for contig_name in group["contig"]:
                if contig_name in fasta_sequences:
                    SeqIO.write(fasta_sequences[contig_name], f_out, "fasta")
        if len(group) == 1:
            single_contig_bins.append((bin_id, group["contig"].iloc[0]))
        else:
            bin_fasta_files.append((bin_id, bin_fasta_file))
    logging.info(f"Extraction terminée : {len(bin_fasta_files)} bins multi-contigs et {len(single_contig_bins)} bins à contig unique.")
    return bin_fasta_files, single_contig_bins, fasta_sequences

def run_flye(bin_id, bin_fasta, flye_dir, flye_threads, logs_dir):
    bin_flye_dir = os.path.join(flye_dir, f"bin_{bin_id}")
    os.makedirs(bin_flye_dir, exist_ok=True)
    log_file = os.path.join(logs_dir, f"flye_bin_{bin_id}.log")

    try:
        with open(log_file, "w") as log:
            subprocess.run([
                "flye",
                "--subassemblies", bin_fasta,
                "--out-dir", bin_flye_dir,
                "--threads", str(flye_threads)
            ], stdout=log, stderr=log, check=True)

        assembled_file = os.path.join(bin_flye_dir, "assembly.fasta")
        if os.path.exists(assembled_file):
            logging.info(f"Assemblage terminé pour le bin {bin_id}.")
            return bin_id, assembled_file
        else:
            raise FileNotFoundError(f"Fichier d'assemblage manquant pour le bin {bin_id}")
    except Exception as e:
        logging.error(f"Échec de l'assemblage pour le bin {bin_id}. Cause probable : {e}")
        return bin_id, bin_fasta  # Renvoie les contigs d'origine en cas d'échec

def run_flye_on_bins(bin_fasta_files, flye_dir, flye_threads, max_parallel, logs_dir):
    logging.info("Début de l'assemblage des bins multi-contigs avec Flye.")
    assembled_files = []
    
    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        future_to_bin = {
            executor.submit(run_flye, bin_id, bin_fasta, flye_dir, flye_threads, logs_dir): (bin_id, bin_fasta)
            for bin_id, bin_fasta in bin_fasta_files
        }
        
        for future in as_completed(future_to_bin):
            bin_id, assembled_file = future.result()
            if assembled_file:
                assembled_files.append((bin_id, assembled_file))
    logging.info(f"Assemblage terminé pour tous les bins multi-contigs.")
    return assembled_files

def process_bin_for_final_file(bin_id, assembly_file, fasta_sequences):
    """
    Prépare les contigs d'un bin pour inclusion dans le fichier final.
    """
    contigs = []
    
    if os.path.exists(assembly_file) and os.path.getsize(assembly_file) > 0:
        contigs.extend(SeqIO.parse(assembly_file, "fasta"))
    else:
        logging.warning(f"Fichier introuvable ou vide pour le bin {bin_id}. Inclusion des contigs d'origine.")
        contigs.extend(SeqIO.parse(assembly_file, "fasta"))
    
    for record in contigs:
        record.id = f"bin_{bin_id}_{record.id}"
        record.description = ""
    
    return contigs


def concatenate_final_contigs(assembled_files, single_contig_bins, fasta_sequences, final_file, threads=4):
    """
    Concatène les contigs finaux (assemblés ou d'origine) dans un fichier unique.
    """
    logging.info("Début de la concaténation des contigs finaux.")
    all_contigs = []
    
    # Gestion en parallèle pour les bins multi-contigs
    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_bin = {
            executor.submit(process_bin_for_final_file, bin_id, assembly_file, fasta_sequences): bin_id
            for bin_id, assembly_file in assembled_files
        }
        
        for future in as_completed(future_to_bin):
            bin_id = future_to_bin[future]
            try:
                contigs = future.result()
                all_contigs.extend(contigs)
            except Exception as e:
                logging.error(f"Erreur lors de la gestion des contigs pour le bin {bin_id} : {e}")
    
    # Ajout des contigs uniques
    logging.info(f"Nombre de bins à contig unique à écrire : {len(single_contig_bins)}")
    for bin_id, contig_id in single_contig_bins:
        if contig_id in fasta_sequences:
            record = fasta_sequences[contig_id]
            record.id = f"bin_{bin_id}_{record.id}"
            record.description = ""
            all_contigs.append(record)
        else:
            logging.warning(f"Contig introuvable : {contig_id} pour le bin {bin_id}.")
    
    # Écriture des contigs dans le fichier final
    with open(final_file, "w") as f_out:
        SeqIO.write(all_contigs, f_out, "fasta")
    
    logging.info(f"Concaténation terminée. Résultat final dans : {final_file}")


def summarize_results(bin_fasta_files, single_contig_bins, assembled_files, final_file):
    """
    Résume les métriques importantes sur les bins et contigs.
    """
    logging.info("Début du résumé des résultats.")

    total_bins = len(bin_fasta_files) + len(single_contig_bins)
    total_initial_contigs = sum(len(list(SeqIO.parse(bin_fasta, "fasta"))) for _, bin_fasta in bin_fasta_files)
    total_initial_contigs += len(single_contig_bins)
    total_final_contigs = sum(1 for _ in SeqIO.parse(final_file, "fasta"))

    # Calcul des réductions de contigs au niveau des bins
    reduced_bins = 0
    total_reduction_bins = 0  # Nombre total de contigs réduits
    for bin_id, assembly_file in assembled_files:
        if os.path.exists(assembly_file) and os.path.getsize(assembly_file) > 0:
            initial_contigs = len(list(SeqIO.parse(assembly_file, "fasta")))
            final_contigs = len(list(SeqIO.parse(assembly_file, "fasta")))
            reduction = initial_contigs - final_contigs
            total_reduction_bins += reduction
            if reduction > 0:
                reduced_bins += 1

    reduction_count = total_initial_contigs - total_final_contigs
    reduction_percent = (reduction_count / total_initial_contigs) * 100 if total_initial_contigs > 0 else 0

    summary_data = {
        "Total bins": total_bins,
        "Multi-contig bins": len(bin_fasta_files),
        "Single-contig bins": len(single_contig_bins),
        "Bins passing Flye": len([f for f in assembled_files if os.path.exists(f[1]) and f[1].endswith(".fasta")]),
        "Failed bins": len([f for f in assembled_files if not os.path.exists(f[1]) or not f[1].endswith(".fasta")]),
        "Reduced bins": reduced_bins,
        "Total reduction in bins (contigs)": total_reduction_bins,
        "Initial contigs": total_initial_contigs,
        "Final contigs": total_final_contigs,
        "Reduction in contigs (count)": reduction_count,
        "Reduction in contigs (%)": f"{reduction_percent:.2f}"
    }

    summary_file = os.path.join(os.path.dirname(final_file), "summary_table.tsv")
    pd.DataFrame(summary_data.items(), columns=["Metric", "Value"]).to_csv(summary_file, sep="\t", index=False)
    logging.info(f"Résumé des résultats sauvegardé dans : {summary_file}")



def main():
    args = parse_arguments()
    bins_dir, flye_dir, final_dir, logs_dir = setup_output_dirs(args.output)
    setup_logging(logs_dir)

    bin_fasta_files, single_contig_bins, fasta_sequences = extract_contigs_by_bin(args.fasta, args.clusters, bins_dir)
    assembled_files = run_flye_on_bins(bin_fasta_files, flye_dir, args.flye_threads, args.max_parallel, logs_dir)

    final_output_file = os.path.join(final_dir, "final_contigs.fasta")
    concatenate_final_contigs(assembled_files, single_contig_bins, fasta_sequences, final_output_file, threads=args.flye_threads)
    summarize_results(bin_fasta_files, single_contig_bins, assembled_files, final_output_file)

    logging.info("Processus terminé.")


if __name__ == "__main__":
    main()
