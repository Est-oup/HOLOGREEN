import os
import sys
import argparse
import logging
from Bio import SeqIO
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Assembler et nettoyer des contigs par bin en supprimant les chevauchements à l'aide de Minimap2.\n\n"
            "STRUCTURE DES ENTRÉES :\n"
            "  1. Fichier FASTA (--fasta) :\n"
            "     - Contient toutes les séquences des contigs.\n"
            "     - Chaque séquence doit avoir un nom unique.\n"
            "\n"
            "  2. Fichier clusters.tsv (--clusters) :\n"
            "     - Fichier TSV (tab-separated values) contenant deux colonnes :\n"
            "         - Colonne 1 : Numéro de bin (identifiant du groupe de contigs).\n"
            "         - Colonne 2 : Nom des contigs (doit correspondre aux noms dans le FASTA).\n"
            "     - Exemple :\n"
            "         1\tcontig_1\n"
            "         1\tcontig_2\n"
            "         2\tcontig_3\n"
            "\n"
            "FONCTIONNALITÉS DU PROGRAMME :\n"
            "  - Regroupe les contigs en bins selon les clusters spécifiés.\n"
            "  - Utilise Minimap2 pour détecter et supprimer les chevauchements entre contigs à l'intérieur de chaque bin.\n"
            "  - Produit des contigs nettoyés pour chaque bin.\n"
            "  - Génère un fichier final concaténant les contigs nettoyés.\n"
            "  - Génère un résumé contenant des statistiques sur les bins et les contigs.\n"
            "\n"
            "SORTIES :\n"
            "  - Un dossier contenant les bins individuels en FASTA.\n"
            "  - Un dossier contenant les contigs nettoyés pour chaque bin.\n"
            "  - Un fichier final : final_contigs.fasta.\n"
            "  - Un résumé des statistiques : summary_table.tsv.\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--fasta", required=True, help="Fichier FASTA contenant les contigs.")
    parser.add_argument("--clusters", required=True, help="Fichier TSV contenant les clusters et contigs (cluster.tsv).")
    parser.add_argument("--output", required=True, help="Dossier de sortie principal.")
    parser.add_argument("--threads", default=1, type=int, help="Nombre de threads pour le traitement (défaut: 1).")
    parser.add_argument("--overlap", default=150, type=int, help="Taille minimale du chevauchement pour supprimer les contigs (défaut: 150).")
    return parser.parse_args()


def setup_output_dirs(output_dir):
    bins_dir = os.path.join(output_dir, "bins")
    results_dir = os.path.join(output_dir, "results")
    final_dir = os.path.join(output_dir, "final_output")
    logs_dir = os.path.join(output_dir, "logs")
    
    for dir_path in [bins_dir, results_dir, final_dir, logs_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    return bins_dir, results_dir, final_dir, logs_dir


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


def run_minimap2(fasta_file, output_paf, overlap_threshold):
    logging.info(f"Exécution de Minimap2 sur {fasta_file} avec seuil de chevauchement {overlap_threshold}...")
    command = f"minimap2 -x ava-pb -m {overlap_threshold} {fasta_file} {fasta_file} > {output_paf}"
    os.system(command)

    if not os.path.exists(output_paf) or os.path.getsize(output_paf) == 0:
        raise RuntimeError(f"Minimap2 a échoué ou le fichier PAF est vide : {output_paf}")
    logging.info(f"Fichier PAF généré avec succès : {output_paf}") 


def parse_paf(paf_file):
    columns = [
        "query_name", "query_length", "query_start", "query_end",
        "strand", "target_name", "target_length", "target_start",
        "target_end", "alignment_length", "mapping_quality"
    ]
    data = []
    with open(paf_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            data.append(fields[:11])  # On garde les 11 premières colonnes
    paf_df = pd.DataFrame(data, columns=columns)

    numeric_cols = ["query_length", "query_start", "query_end", 
                    "target_length", "target_start", "target_end", "alignment_length"]
    paf_df[numeric_cols] = paf_df[numeric_cols].astype(int)

    paf_df = paf_df[paf_df["alignment_length"] >= 150]
    return paf_df


def filter_overlapping_contigs(fasta_file, paf_df, output_file, overlap_threshold):
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    contigs_to_keep = set(sequences.keys())
    modified_sequences = {}

    for _, row in paf_df.iterrows():
        query, target = row["query_name"], row["target_name"]
        if query == target:
            continue

        query_overlap_end = row["query_start"] <= overlap_threshold
        target_overlap_start = row["target_end"] >= row["target_length"] - overlap_threshold

        if query_overlap_end and target_overlap_start:
            logging.info(f"Chevauchement détecté entre {query} et {target}.")
            if row["query_length"] >= row["target_length"]:
                contigs_to_keep.discard(target)
            else:
                contigs_to_keep.discard(query)

    with open(output_file, "w") as f_out:
        for contig_id in contigs_to_keep:
            f_out.write(f">{contig_id}\n{sequences[contig_id]}\n")

    logging.info(f"Fichier FASTA ajusté sauvegardé dans : {output_file}")


def process_bin(bin_id, bin_fasta, results_dir, overlap_threshold):
    paf_file = os.path.join(results_dir, f"bin_{bin_id}.paf")
    filtered_fasta = os.path.join(results_dir, f"filtered_bin_{bin_id}.fasta")

    try:
        run_minimap2(bin_fasta, paf_file, overlap_threshold)
        paf_df = parse_paf(paf_file)
        filter_overlapping_contigs(bin_fasta, paf_df, filtered_fasta, overlap_threshold)
        return bin_id, filtered_fasta, True
    except Exception as e:
        logging.error(f"Erreur lors du traitement du bin {bin_id} : {e}")
        return bin_id, bin_fasta, False


def process_bins_in_parallel(bin_fasta_files, results_dir, overlap_threshold, threads):
    logging.info("Début du traitement des bins pour supprimer les chevauchements.")
    filtered_bins = []
    successful_bins = 0
    failed_bins = 0

    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_bin = {
            executor.submit(process_bin, bin_id, bin_fasta, results_dir, overlap_threshold): (bin_id, bin_fasta)
            for bin_id, bin_fasta in bin_fasta_files
        }
        
        for future in as_completed(future_to_bin):
            bin_id, filtered_fasta, success = future.result()
            if success:
                successful_bins += 1
                logging.info(f"Bin {bin_id} traité avec succès.")
            else:
                failed_bins += 1
                logging.warning(f"Bin {bin_id} a échoué. Utilisation des contigs d'origine.")
            filtered_bins.append((bin_id, filtered_fasta))
    
    logging.info(f"Traitement des bins terminé : {successful_bins} réussis, {failed_bins} échoués.")
    return filtered_bins, successful_bins, failed_bins


def concatenate_final_contigs(filtered_bins, single_contig_bins, fasta_sequences, final_file):
    logging.info("Concaténation des contigs finaux.")
    all_contigs = []

    for bin_id, filtered_fasta in filtered_bins:
        all_contigs.extend(SeqIO.parse(filtered_fasta, "fasta"))

    for bin_id, contig_id in single_contig_bins:
        if contig_id in fasta_sequences:
            record = fasta_sequences[contig_id]
            record.id = f"bin_{bin_id}_{record.id}"
            record.description = ""
            all_contigs.append(record)

    with open(final_file, "w") as f_out:
        SeqIO.write(all_contigs, f_out, "fasta")

    logging.info(f"Fichier final écrit dans : {final_file}")


def summarize_results(bin_fasta_files, single_contig_bins, filtered_bins, final_file, summary_file, successful_bins, failed_bins):
    total_bins = len(bin_fasta_files) + len(single_contig_bins)
    total_initial_contigs = sum(len(list(SeqIO.parse(bin_fasta, "fasta"))) for _, bin_fasta in bin_fasta_files)
    total_initial_contigs += len(single_contig_bins)
    total_final_contigs = sum(1 for _ in SeqIO.parse(final_file, "fasta"))
    total_initial_bases = sum(len(record.seq) for record in SeqIO.parse(final_file, "fasta"))

    reduction_count = total_initial_contigs - total_final_contigs
    reduction_percent = (reduction_count / total_initial_contigs) * 100 if total_initial_contigs > 0 else 0

    summary_data = {
        "Total bins": total_bins,
        "Successful bins": successful_bins,
        "Failed bins": failed_bins,
        "Initial contigs": total_initial_contigs,
        "Final contigs": total_final_contigs,
        "Reduction in contigs (count)": reduction_count,
        "Reduction in contigs (%)": f"{reduction_percent:.2f}",
        "Total bases in final contigs": total_initial_bases
    }

    pd.DataFrame(summary_data.items(), columns=["Metric", "Value"]).to_csv(summary_file, sep="\t", index=False)
    logging.info(f"Résumé des résultats sauvegardé dans : {summary_file}")


def main():
    args = parse_arguments()
    bins_dir, results_dir, final_dir, logs_dir = setup_output_dirs(args.output)
    setup_logging(logs_dir)

    try:
        # Extraction des contigs par bin
        bin_fasta_files, single_contig_bins, fasta_sequences = extract_contigs_by_bin(args.fasta, args.clusters, bins_dir)

        # Traitement des bins en parallèle
        filtered_bins, successful_bins, failed_bins = process_bins_in_parallel(
            bin_fasta_files, results_dir, args.overlap, args.threads
        )
        
        # Concaténation des résultats finaux
        final_output_file = os.path.join(final_dir, "final_contigs.fasta")
        concatenate_final_contigs(filtered_bins, single_contig_bins, fasta_sequences, final_output_file)

        # Génération du résumé des résultats
        summary_file = os.path.join(final_dir, "summary_table.tsv")
        summarize_results(
            bin_fasta_files, single_contig_bins, filtered_bins, final_output_file, summary_file, successful_bins, failed_bins
        )

        logging.info("Processus terminé avec succès.")
    except Exception as e:
        logging.error(f"Échec global du processus : {e}")
        logging.info("Les contigs d'origine seront utilisés.")

        # Fichier final avec les contigs d'origine
        final_output_file = os.path.join(final_dir, "original_contigs.fasta")
        with open(final_output_file, "w") as f_out:
            SeqIO.write(SeqIO.parse(args.fasta, "fasta"), f_out, "fasta")

        # Résumé minimal en cas d'échec
        summary_file = os.path.join(final_dir, "summary_table.tsv")
        summary_data = {
            "Total bins": "N/A",
            "Successful bins": "N/A",
            "Failed bins": "N/A",
            "Initial contigs": len(list(SeqIO.parse(args.fasta, "fasta"))),
            "Final contigs": "N/A",
            "Reduction in contigs (count)": "N/A",
            "Reduction in contigs (%)": "N/A",
            "Total bases in final contigs": "N/A"
        }
        pd.DataFrame(summary_data.items(), columns=["Metric", "Value"]).to_csv(summary_file, sep="\t", index=False)

        logging.info(f"Résumé minimal des résultats sauvegardé dans : {summary_file}")
        logging.info(f"Fichier final écrit avec les contigs originaux : {final_output_file}")


if __name__ == "__main__":
    main()

