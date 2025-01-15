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
        description="Assembler des contigs par bin avec Minimap2 et Racon."
    )
    parser.add_argument("--fasta", required=True, help="Fichier FASTA contenant les contigs.")
    parser.add_argument("--clusters", required=True, help="Fichier TSV contenant les clusters et contigs (cluster.tsv).")
    parser.add_argument("--output", required=True, help="Dossier de sortie principal.")
    parser.add_argument("--threads", default=4, type=int, help="Nombre de threads pour Minimap2 et Racon (défaut: 4).")
    parser.add_argument("--max_parallel", default=4, type=int, help="Nombre maximum de bins à assembler en parallèle (défaut: 4).")
    parser.add_argument("--minimap_x", default="asm5", type=str, help="Paramètres généraux de Minimap2 (défaut: asm5).")
    
    return parser.parse_args()


def setup_output_dirs(output_dir):
    bins_dir = os.path.join(output_dir, "bins")
    racon_dir = os.path.join(output_dir, "racon_results")
    final_dir = os.path.join(output_dir, "final_output")
    logs_dir = os.path.join(output_dir, "logs")
    
    for dir_path in [bins_dir, racon_dir, final_dir, logs_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    return bins_dir, racon_dir, final_dir, logs_dir


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

    multi_contig_bins = []
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
            multi_contig_bins.append((bin_id, bin_fasta_file))
    logging.info(f"Extraction terminée : {len(multi_contig_bins)} bins multi-contigs et {len(single_contig_bins)} bins à contig unique.")
    return multi_contig_bins, single_contig_bins, fasta_sequences


def run_minimap2_racon(bin_id, bin_fasta, racon_dir, threads, logs_dir, minimap_x):
    bin_racon_dir = os.path.join(racon_dir, f"bin_{bin_id}")
    os.makedirs(bin_racon_dir, exist_ok=True)
    log_file = os.path.join(logs_dir, f"minimap2_racon_bin_{bin_id}.log")
    
    paf_file = os.path.join(bin_racon_dir, "overlaps.paf")
    racon_output_file = os.path.join(bin_racon_dir, "racon_out.fasta")

    try:
        # Étape 1 : Minimap2
        with open(log_file, "w") as log:
            result = subprocess.run(
                ["minimap2", "-x", str(minimap_x), "-t", str(threads), bin_fasta, bin_fasta],
                stdout=open(paf_file, "w"), stderr=log
            )
            if result.returncode != 0 or os.path.getsize(paf_file) == 0:
                logging.warning(f"Minimap2 a échoué pour le bin {bin_id}.")
                return bin_id, bin_fasta  # Retourne le fichier d'origine

        # Étape 2 : Racon
        with open(log_file, "a") as log:
            result = subprocess.run(
                ["racon", "--no-trimming", "-t", str(threads), bin_fasta, paf_file, bin_fasta],
                stdout=open(racon_output_file, "w"), stderr=log
            )
            if result.returncode != 0 or os.path.getsize(racon_output_file) == 0:
                logging.warning(f"Racon a échoué pour le bin {bin_id}.")
                return bin_id, bin_fasta  # Retourne le fichier d'origine

        logging.info(f"Assemblage réussi pour le bin {bin_id}.")
        return bin_id, racon_output_file

    except Exception as e:
        logging.error(f"Échec pour le bin {bin_id}. Erreur : {e}")
        return bin_id, bin_fasta  # Retourne le fichier d'origine


def run_minimap2_racon_on_bins(multi_contig_bins, racon_dir, threads, max_parallel, logs_dir, minimap_x):
    logging.info("Début de l'assemblage des bins multi-contigs avec Minimap2 et Racon.")
    results = []

    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        future_to_bin = {
            executor.submit(run_minimap2_racon, bin_id, bin_fasta, racon_dir, threads, logs_dir, minimap_x): bin_id
            for bin_id, bin_fasta in multi_contig_bins
        }
        
        for future in as_completed(future_to_bin):
            bin_id, result = future.result()
            results.append((bin_id, result))
    logging.info("Assemblage terminé.")
    return results


def process_bin_for_final_file(bin_id, assembly_file, fasta_sequences):
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


def concatenate_final_contigs(multi_contig_results, single_contig_bins, fasta_sequences, final_file, threads=4):
    logging.info("Début de la concaténation des contigs finaux.")
    all_contigs = []
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_bin = {
            executor.submit(process_bin_for_final_file, bin_id, assembly_file, fasta_sequences): bin_id
            for bin_id, assembly_file in multi_contig_results
        }
        for future in as_completed(future_to_bin):
            try:
                contigs = future.result()
                all_contigs.extend(contigs)
            except Exception as e:
                logging.error(f"Erreur lors de la gestion des contigs : {e}")
    
    for bin_id, contig_id in single_contig_bins:
        if contig_id in fasta_sequences:
            record = fasta_sequences[contig_id]
            record.id = f"bin_{bin_id}_{record.id}"
            record.description = ""
            all_contigs.append(record)
        else:
            logging.warning(f"Contig introuvable : {contig_id}")
    
    with open(final_file, "w") as f_out:
        SeqIO.write(all_contigs, f_out, "fasta")
    logging.info(f"Concaténation terminée. Résultat dans : {final_file}")


def summarize_results(multi_contig_bins, single_contig_bins, multi_contig_results, final_file):
    logging.info("Début du résumé des résultats.")

    # Nombre total de bins
    total_bins = len(multi_contig_bins) + len(single_contig_bins)

    # Nombre total de contigs initiaux
    total_contigs_initial = sum(len(list(SeqIO.parse(bin_fasta, "fasta"))) for _, bin_fasta in multi_contig_bins)
    total_contigs_initial += len(single_contig_bins)

    # Nombre total de contigs finaux
    total_contigs_final = sum(1 for _ in SeqIO.parse(final_file, "fasta"))

    # Comptage des bins échoués et réussis
    minimap_failures = len([result for _, result in multi_contig_results if result == "Minimap_failed"])
    racon_failures = len([result for _, result in multi_contig_results if result == "Racon_failed"])
    successful_bins = len([result for _, result in multi_contig_results if result not in ["Minimap_failed", "Racon_failed"]])

    # Réduction de contigs
    reduction_count = total_contigs_initial - total_contigs_final
    reduction_percent = (reduction_count / total_contigs_initial * 100) if total_contigs_initial > 0 else 0

    # Création d'un tableau de résumé
    summary_data = {
        "Total bins": total_bins,
        "Multi-contig bins": len(multi_contig_bins),
        "Single-contig bins": len(single_contig_bins),
        "Bins passing Minimap2 + Racon": successful_bins,
        "Bins échoués (Minimap2)": minimap_failures,
        "Bins échoués (Racon)": racon_failures,
        "Total initial contigs": total_contigs_initial,
        "Final contig count": total_contigs_final,
        "Reduction in contigs (count)": reduction_count,
        "Reduction in contigs (%)": f"{reduction_percent:.2f}"
    }

    # Sauvegarde du tableau dans un fichier TSV
    summary_file = os.path.join(os.path.dirname(final_file), "summary_table.tsv")
    pd.DataFrame(summary_data.items(), columns=["Metric", "Value"]).to_csv(summary_file, sep="\t", index=False)

    logging.info(f"Résumé des résultats sauvegardé dans : {summary_file}")



def main():
    args = parse_arguments()
    bins_dir, racon_dir, final_dir, logs_dir = setup_output_dirs(args.output)
    setup_logging(logs_dir)

    multi_contig_bins, single_contig_bins, fasta_sequences = extract_contigs_by_bin(args.fasta, args.clusters, bins_dir)
    multi_contig_results = run_minimap2_racon_on_bins(multi_contig_bins, racon_dir, args.threads, args.max_parallel, logs_dir, args.minimap_x)

    final_output_file = os.path.join(final_dir, "final_contigs.fasta")
    concatenate_final_contigs(multi_contig_results, single_contig_bins, fasta_sequences, final_output_file, threads=args.threads)
    summarize_results(multi_contig_bins, single_contig_bins, multi_contig_results, final_output_file)

    logging.info("Processus terminé.")


if __name__ == "__main__":
    main()
