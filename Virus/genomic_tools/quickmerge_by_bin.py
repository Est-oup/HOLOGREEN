import os
import sys
import argparse
from Bio import SeqIO
import pandas as pd
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configuration des logs
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("process.log"),
        logging.StreamHandler(sys.stdout)
    ]
)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Assembler des contigs par bin avec QuickMerge."
    )
    parser.add_argument("--fasta", required=True, help="Fichier FASTA contenant les contigs.")
    parser.add_argument("--clusters", required=True, help="Fichier TSV contenant les clusters et contigs (cluster.tsv).")
    parser.add_argument("--output", required=True, help="Dossier de sortie principal.")
    parser.add_argument("--threads", default=4, type=int, help="Nombre de threads pour traiter les bins en parallèle (défaut: 4).")
    parser.add_argument("--hco", default="5.0", help="High Confidence Overlap pour QuickMerge (défaut: 5.0).")
    parser.add_argument("--c", default="1.5", help="Minimum contig length ratio pour QuickMerge (défaut: 1.5).")
    parser.add_argument("--seed_length", default="1000", help="Longueur minimale des contigs pour QuickMerge (défaut: 1000).")
    parser.add_argument("--merge_length", default="5000", help="Longueur minimale des overlaps pour QuickMerge (défaut: 5000).")
    return parser.parse_args()

def setup_output_dirs(output_dir):
    res_bin_dir = os.path.join(output_dir, "res_bin")
    final_fasta_dir = os.path.join(output_dir, "final_fasta")
    os.makedirs(res_bin_dir, exist_ok=True)
    os.makedirs(final_fasta_dir, exist_ok=True)
    return res_bin_dir, final_fasta_dir

def extract_contigs_by_bin(fasta_file, cluster_file, res_bin_dir):
    logging.info("Début de l'extraction des contigs par bin.")
    clusters = pd.read_csv(cluster_file, sep="\t", header=None, names=["bin", "contig"])
    fasta_sequences = {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}
    bins = clusters.groupby("bin")

    bin_fasta_files = []
    single_contig_bins = []
    for bin_id, group in bins:
        bin_dir = os.path.join(res_bin_dir, f"bin_{bin_id}")
        os.makedirs(bin_dir, exist_ok=True)
        
        bin_fasta_file = os.path.join(bin_dir, f"bin_{bin_id}_initial.fasta")
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

def run_quickmerge(bin_id, bin_fasta, res_bin_dir, hco, c, seed_length, merge_length):
    bin_dir = os.path.join(res_bin_dir, f"bin_{bin_id}")
    log_file = os.path.join(bin_dir, "quickmerge.log")
    delta_file = os.path.join(bin_dir, "temp_delta.delta")
    filtered_delta_file = os.path.join(bin_dir, "filtered_delta.delta")
    quickmerge_output_prefix = os.path.join(bin_dir, f"bin_{bin_id}_merged_")
    quickmerge_output_file = f"{quickmerge_output_prefix}fasta"

    try:
        # Étape 1 : Générer le fichier delta avec MUMmer
        subprocess.run(
            ["nucmer", "--maxmatch", "-l", "200", "-prefix", os.path.join(bin_dir, "temp_delta"), bin_fasta, bin_fasta],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=True
        )

        # Étape 2 : Filtrer les alignements avec delta-filter
        subprocess.run(
            ["delta-filter", "-l", "200", delta_file],
            stdout=open(filtered_delta_file, "w"), stderr=subprocess.PIPE, check=True
        )

        # Étape 3 : Lancer QuickMerge
        result = subprocess.run(
            [
                "quickmerge",
                "-d", filtered_delta_file,
                "-q", bin_fasta,
                "-r", bin_fasta,
                "-hco", hco,
                "-c", c,
                "-l", seed_length,
                "-ml", merge_length,
                "-p", quickmerge_output_prefix,
            ],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
        )

        # Afficher uniquement les résultats pertinents si succès
        if os.path.exists(quickmerge_output_file) and os.path.getsize(quickmerge_output_file) > 0:
            logging.info(f"Fusion réussie pour le bin {bin_id}.")
            return bin_id, True
        else:
            logging.warning(f"Aucun résultat généré pour le bin {bin_id}.")
            return bin_id, False

    except subprocess.CalledProcessError as e:
        logging.error(f"QuickMerge a échoué pour le bin {bin_id}. Erreur : {e.stderr.strip()}")
        return bin_id, False




def run_quickmerge_on_bins(bin_fasta_files, res_bin_dir, hco, c, seed_length, merge_length, max_parallel):
    logging.info("Début de l'assemblage des bins multi-contigs avec QuickMerge.")
    assembled_files = []

    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        future_to_bin = {
            executor.submit(run_quickmerge, bin_id, bin_fasta, res_bin_dir, hco, c, seed_length, merge_length): (bin_id, bin_fasta)
            for bin_id, bin_fasta in bin_fasta_files
        }
        
        for future in as_completed(future_to_bin):
            bin_id, assembled_file = future.result()
            if assembled_file:
                assembled_files.append((bin_id, assembled_file))
    logging.info(f"Assemblage terminé pour tous les bins multi-contigs.")
    return assembled_files

def concatenate_final_contigs(assembled_files, single_contig_bins, fasta_sequences, final_file):
    logging.info("Début de la concaténation des contigs finaux.")
    with open(final_file, "w") as f_out:
        for bin_id, assembly_file in assembled_files:
            for record in SeqIO.parse(assembly_file, "fasta"):
                record.id = f"bin_{bin_id}_{record.id}"
                record.description = ""
                SeqIO.write(record, f_out, "fasta")
        for bin_id, contig_id in single_contig_bins:
            record = fasta_sequences[contig_id]
            record.id = f"bin_{bin_id}_{record.id}"
            record.description = ""
            SeqIO.write(record, f_out, "fasta")
    logging.info(f"Concaténation terminée. Résultat final dans : {final_file}")

def main():
    args = parse_arguments()
    res_bin_dir, final_fasta_dir = setup_output_dirs(args.output)
    bin_fasta_files, single_contig_bins, fasta_sequences = extract_contigs_by_bin(args.fasta, args.clusters, res_bin_dir)
    assembled_files = run_quickmerge_on_bins(
        bin_fasta_files, res_bin_dir, args.hco, args.c, args.seed_length, args.merge_length, args.threads
    )
    final_output_file = os.path.join(final_fasta_dir, "final_contigs.fasta")
    concatenate_final_contigs(assembled_files, single_contig_bins, fasta_sequences, final_output_file)
    logging.info("Processus terminé.")

if __name__ == "__main__":
    main()
