import os
import sys
import pandas as pd
from pathlib import Path

def read_worst_scores(worst_scores_folder):
    worst_scores = {}
    for file in Path(worst_scores_folder).glob("*_worst_scores.txt"):
        df = pd.read_csv(file, sep="\t", header=None)
        # Supposons que la dernière colonne contient les scores
        worst_score = df.iloc[:, -1].min()
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
        gene_name = file.stem.replace("_filtered", "").replace("-polinto", "")
        if gene_name not in worst_scores:
            continue
        score_threshold = worst_scores[gene_name]
        df = pd.read_csv(file, sep="\t", header=None)
        
        # Filtrer par score
        filtered_df = df[df.iloc[:, -2] <= score_threshold]  # Avant-dernière colonne pour les scores
        output_file = Path(output_folder) / file.name
        filtered_df.to_csv(output_file, sep="\t", index=False, header=False)

def extract_contig_ids(alignment_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    for file in Path(alignment_folder).glob("*"):
        df = pd.read_csv(file, sep="\t", header=None)
        
        if "filtered" in file.name:
            contig_ids = df.iloc[:, 0].unique()  # Première colonne pour les contigs dans les fichiers filtered
        else:
            contig_ids = df.iloc[:, 1].unique()  # Deuxième colonne pour les contigs dans les fichiers non filtered
        
        output_file = Path(output_folder) / f"{file.stem}.txt"
        pd.Series(contig_ids).to_csv(output_file, index=False, header=False)

def extract_fasta_sequences(fasta_folder, contig_ids_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    for contig_file in Path(contig_ids_folder).glob("*.txt"):
        gene_name = contig_file.stem.replace("", "")
        fasta_file = Path(fasta_folder) / f"{gene_name}.fasta"
        
        if not fasta_file.exists():
            continue
        
        # Lire les IDs des contigs
        contig_ids = set(pd.read_csv(contig_file, header=None).iloc[:, 0].tolist())
        
        # Lire le fichier FASTA et écrire les séquences correspondant aux contigs
        with open(fasta_file, "r") as fasta, open(Path(output_folder) / f"{gene_name}.fasta", "w") as output_fasta:
            write = False
            for line in fasta:
                if line.startswith(">"):
                    contig_id = line[1:].strip().split()[0]  # Extraire l'ID de contig
                    write = contig_id in contig_ids
                if write:
                    output_fasta.write(line)

def main():
    worst_scores_folder = sys.argv[1]
    alignment_folder = sys.argv[2]
    output_alignment_folder = sys.argv[3]
    contig_ids_folder = sys.argv[4]
    fasta_folder = sys.argv[5]
    output_fasta_folder = sys.argv[6]
    summary_file = sys.argv[7]

    worst_scores = read_worst_scores(worst_scores_folder)
    write_summary(worst_scores, summary_file)
    filter_alignments_by_score(alignment_folder, worst_scores, output_alignment_folder)
    extract_contig_ids(output_alignment_folder, contig_ids_folder)
    extract_fasta_sequences(fasta_folder, contig_ids_folder, output_fasta_folder)

if __name__ == "__main__":
    main()
