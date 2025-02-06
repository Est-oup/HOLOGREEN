import os
import sys
import pandas as pd
from pathlib import Path

def read_worst_scores(worst_scores_folder):
    worst_scores = {}
    for file in Path(worst_scores_folder).glob("*_worst_scores.txt"):
        df = pd.read_csv(file, sep="\t", header=None)
        worst_score = df.iloc[:, -1].min()  # Prendre le score minimum
        gene_name = file.stem.replace("_worst_scores", "")
        worst_scores[gene_name] = worst_score
    return worst_scores

def write_summary(worst_scores, summary_file):
    with open(summary_file, "w") as f:
        f.write("Gene_Marker\tWorst_Score\n")
        for gene, score in worst_scores.items():
            f.write(f"{gene}\t{score}\n")

def filter_alignments_by_score(alignment_folder, worst_scores, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    for file in Path(alignment_folder).glob("*"):
        gene_name = file.stem.split("_")[0]  # Extraire uniquement le nom du marqueur avant le premier underscore

        if gene_name not in worst_scores:
            print(f"Skipping {file.name}, no worst score available.")
            continue

        score_threshold = worst_scores[gene_name]
        df = pd.read_csv(file, sep="\t", header=None)

        # Filtrer les alignements par score
        filtered_df = df[df.iloc[:, -2] <= score_threshold]  # Avant-dernière colonne pour les scores
        print(f"{len(filtered_df)} alignments retained for {gene_name} (threshold: {score_threshold})")

        output_file = Path(output_folder) / file.name
        filtered_df.to_csv(output_file, sep="\t", index=False, header=False)

def extract_contig_ids(alignment_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    for file in Path(alignment_folder).glob("*"):
        df = pd.read_csv(file, sep="\t", header=None)
        contig_ids = df.iloc[:, 0].unique()  # Première colonne pour les contigs dans les fichiers filtrés
        output_file = Path(output_folder) / f"{file.stem}.txt"
        pd.Series(contig_ids).to_csv(output_file, index=False, header=False)
        print(f"Extracted {len(contig_ids)} contig IDs for {file.name}")

def find_fasta_file(fasta_folder, gene_name):
    """ Trouver le fichier FASTA correspondant en vérifiant d'abord sans _filtered puis avec. """
    non_filtered_file = Path(fasta_folder) / f"{gene_name}.fasta"
    filtered_file = Path(fasta_folder) / f"{gene_name}_filtered.fasta"
    
    if non_filtered_file.exists():
        print(f"Using non-filtered file: {non_filtered_file}")
        return non_filtered_file
    elif filtered_file.exists():
        print(f"Using filtered file: {filtered_file}")
        return filtered_file
    else:
        print(f"FASTA file not found for {gene_name}")
        return None

def extract_fasta_sequences(fasta_folder, contig_ids_folder, output_folder):
    """Extraire les séquences FASTA en fonction des identifiants de contigs."""
    os.makedirs(output_folder, exist_ok=True)
    
    for contig_file in Path(contig_ids_folder).glob("*.txt"):
        # Utiliser le nom complet du fichier sans modification
        gene_name = contig_file.stem
        
        # Trouver le fichier FASTA correspondant avec le même nom de base
        fasta_file = find_fasta_file(fasta_folder, gene_name)

        if not fasta_file:
            continue

        # Lire les IDs des contigs
        contig_ids = set(pd.read_csv(contig_file, header=None).iloc[:, 0].tolist())

        # Lire le fichier FASTA et extraire les séquences correspondant aux contigs
        found_contigs = 0
        with open(fasta_file, "r") as fasta, open(Path(output_folder) / f"{gene_name}.fasta", "w") as output_fasta:
            write = False
            for line in fasta:
                if line.startswith(">"):
                    contig_id = line[1:].split()[0].strip()  # Extraire l'ID de contig
                    write = contig_id in contig_ids  # Vérifier si l'ID de contig est dans la liste des contigs extraits
                    if write:
                        found_contigs += 1
                if write:
                    output_fasta.write(line)

        # Log des contigs extraits
        print(f"Extracted {found_contigs} contigs for {gene_name} from {fasta_file}")

def main():
    worst_scores_folder = sys.argv[1]
    alignment_folder = sys.argv[2]
    output_alignment_folder = sys.argv[3]
    contig_ids_folder = sys.argv[4]
    fasta_folder = sys.argv[5]
    output_fasta_folder = sys.argv[6]
    summary_file = sys.argv[7]

    print(f"Reading worst scores from {worst_scores_folder}")
    worst_scores = read_worst_scores(worst_scores_folder)
    
    print(f"Writing summary of worst scores to {summary_file}")
    write_summary(worst_scores, summary_file)
    
    print(f"Filtering alignments by score from {alignment_folder}")
    filter_alignments_by_score(alignment_folder, worst_scores, output_alignment_folder)
    
    print(f"Extracting contig IDs to {contig_ids_folder}")
    extract_contig_ids(output_alignment_folder, contig_ids_folder)
    
    print(f"Extracting FASTA sequences based on contig IDs from {contig_ids_folder} and saving to {output_fasta_folder}")
    extract_fasta_sequences(fasta_folder, contig_ids_folder, output_fasta_folder)

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python STEP2C_filter_by_ref_score.py <worst_scores_folder> <alignment_folder> <output_alignment_folder> <contig_ids_folder> <fasta_folder> <output_fasta_folder> <summary_file>")
        sys.exit(1)

    main()
