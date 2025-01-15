import os
import pandas as pd
from Bio import SeqIO

def run_minimap2(fasta_file, output_paf, overlap_threshold):
    """
    Exécute Minimap2 pour détecter les chevauchements entre contigs.
    """
    print("Exécution de Minimap2...")
    command = f"minimap2 -x ava-ont -m {overlap_threshold} {fasta_file} {fasta_file} > {output_paf}"
    os.system(command)
    print(f"Fichier PAF généré : {output_paf}")

def parse_paf(paf_file):
    """
    Parse le fichier PAF en DataFrame.
    """
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

    # Convertir les colonnes nécessaires en entiers
    numeric_cols = ["query_length", "query_start", "query_end", 
                    "target_length", "target_start", "target_end", "alignment_length"]
    paf_df[numeric_cols] = paf_df[numeric_cols].astype(int)

    # Filtrer les alignements courts
    paf_df = paf_df[paf_df["alignment_length"] >= 150]

    return paf_df

def filter_overlapping_contigs(fasta_file, paf_df, output_file, overlap_threshold=150):
    """
    Supprime un des deux contigs impliqués dans un chevauchement.
    """
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    contigs_to_keep = set(sequences.keys())

    for _, row in paf_df.iterrows():
        query, target = row["query_name"], row["target_name"]
        if query == target:
            continue  # Ignorer les alignements internes

        # Vérifier si les chevauchements sont aux extrémités
        query_overlap_end = row["query_start"] <= overlap_threshold
        target_overlap_start = row["target_end"] >= row["target_length"] - overlap_threshold

        if query_overlap_end and target_overlap_start:
            print(f"Chevauchement détecté : {query} avec {target}.")
            # Supprimer le contig le plus court
            if row["query_length"] >= row["target_length"]:
                contigs_to_keep.discard(target)
            else:
                contigs_to_keep.discard(query)

    # Sauvegarder les contigs restants dans un fichier FASTA
    with open(output_file, "w") as f_out:
        for contig_id in contigs_to_keep:
            f_out.write(f">{contig_id}\n{sequences[contig_id]}\n")

    print(f"Fichier FASTA final sauvegardé dans : {output_file}")

def main():
    fasta_file = "bin_71.fasta"  # Fichier FASTA d'entrée
    output_paf = "overlaps.paf"  # Fichier PAF généré par Minimap2
    output_file = "final_contigs.fasta"  # Fichier FASTA final
    overlap_threshold = 100  # Seuil minimal de chevauchement

    # Étape 1 : Exécuter Minimap2 pour générer le fichier PAF
    run_minimap2(fasta_file, output_paf, overlap_threshold)

    # Étape 2 : Analyser le PAF et supprimer les contigs en double
    paf_df = parse_paf(output_paf)
    filter_overlapping_contigs(fasta_file, paf_df, output_file, overlap_threshold)

if __name__ == "__main__":
    main()
