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
        description="Assembler des contigs par bin avec Flye et organiser les sorties."
    )
    parser.add_argument("--fasta", required=True, help="Fichier FASTA contenant les contigs.")
    parser.add_argument("--clusters", required=True, help="Fichier TSV contenant les clusters et contigs (cluster.tsv).")
    parser.add_argument("--output", required=True, help="Dossier de sortie principal.")
    parser.add_argument("--flye_threads", default=4, type=int, help="Nombre de threads pour Flye (défaut: 4).")
    parser.add_argument("--max_parallel", default=4, type=int, help="Nombre maximum de bins à assembler en parallèle (défaut: 4).")
    return parser.parse_args()

def setup_output_dirs(output_dir):
    bins_dir = os.path.join(output_dir, "bins")
    flye_dir = os.path.join(output_dir, "flye_results")
    final_dir = os.path.join(output_dir, "final_output")
    logs_dir = os.path.join(output_dir, "logs")
    benchmark_file = os.path.join(output_dir, "benchmark.csv")
    
    for dir_path in [bins_dir, flye_dir, final_dir, logs_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    return bins_dir, flye_dir, final_dir, logs_dir, benchmark_file

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
    return bin_fasta_files, single_contig_bins, clusters

def calculate_genome_size(bin_fasta_file):
    total_size = 0
    for record in SeqIO.parse(bin_fasta_file, "fasta"):
        total_size += len(record.seq)
    return max(int(total_size * 1.2), 100000)

def run_flye(bin_id, bin_fasta, flye_dir, flye_threads, logs_dir):
    bin_flye_dir = os.path.join(flye_dir, f"bin_{bin_id}")
    os.makedirs(bin_flye_dir, exist_ok=True)
    log_file = os.path.join(logs_dir, f"flye_bin_{bin_id}.log")
    genome_size = calculate_genome_size(bin_fasta)

    try:
        with open(log_file, "w") as log:
            subprocess.run([
                "flye",
                "--subassemblies", bin_fasta,
                "--genome-size", f"{genome_size}",
                "--out-dir", bin_flye_dir,
                "--threads", str(flye_threads)
            ], stdout=log, stderr=log, check=True)

        assembled_file = os.path.join(bin_flye_dir, "assembly.fasta")
        if os.path.exists(assembled_file):
            logging.info(f"Assemblage terminé pour le bin {bin_id} (taille estimée : {genome_size}).")
            return bin_id, assembled_file
        else:
            raise FileNotFoundError(f"Fichier d'assemblage manquant pour le bin {bin_id}")
    except Exception as e:
        logging.error(f"Échec de l'assemblage pour le bin {bin_id}. Contigs d'origine utilisés.")
        return bin_id, bin_fasta

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

def benchmark_results(clusters, final_file, benchmark_file):
    logging.info("Création du fichier de benchmark.")
    bins = clusters.groupby("bin")
    final_contigs = [record.id for record in SeqIO.parse(final_file, "fasta")]

    data = []
    for bin_id, group in bins:
        original_count = len(group)
        final_count = len([contig for contig in final_contigs if contig.startswith(f"bin_{bin_id}_")])
        data.append({
            "bin": bin_id,
            "contigs_initial": original_count,
            "contigs_final": final_count,
            "difference": final_count - original_count
        })

    df = pd.DataFrame(data)
    df.to_csv(benchmark_file, index=False)
    logging.info(f"Fichier de benchmark créé : {benchmark_file}")

def main():
    args = parse_arguments()
    bins_dir, flye_dir, final_dir, logs_dir, benchmark_file = setup_output_dirs(args.output)
    
    bin_fasta_files, single_contig_bins, clusters = extract_contigs_by_bin(args.fasta, args.clusters, bins_dir)
    fasta_sequences = {record.id: record for record in SeqIO.parse(args.fasta, "fasta")}
    
    assembled_files = run_flye_on_bins(bin_fasta_files, flye_dir, args.flye_threads, args.max_parallel, logs_dir)
    
    final_output_file = os.path.join(final_dir, "final_contigs.fasta")
    concatenate_final_contigs(assembled_files, single_contig_bins, fasta_sequences, final_output_file)
    
    benchmark_results(clusters, final_output_file, benchmark_file)
    
    logging.info("Processus terminé.")

if __name__ == "__main__":
    main()
