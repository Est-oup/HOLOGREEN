import sys
from Bio import SeqIO
import pandas as pd

def extract_contigs(bin_tsv, fasta_file, target_bin, output_fasta):
    # Charger le fichier TSV avec les bins et contigs
    bin_data = pd.read_csv(bin_tsv, sep="\t", header=None, names=["bin", "contig"])

    # Filtrer les contigs correspondant au bin cible
    contigs_to_extract = bin_data[bin_data["bin"] == int(target_bin)]["contig"]

    if contigs_to_extract.empty:
        print(f"Aucun contig trouvé pour le bin {target_bin}.")
        return

    # Charger les séquences du fichier FASTA
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # Écrire les contigs correspondants dans le fichier de sortie
    with open(output_fasta, "w") as f_out:
        for contig in contigs_to_extract:
            if contig in fasta_sequences:
                SeqIO.write(fasta_sequences[contig], f_out, "fasta")
            else:
                print(f"Contig {contig} non trouvé dans le fichier FASTA.")
    print(f"Extraction terminée pour le bin {target_bin}. Résultats enregistrés dans {output_fasta}.")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python extract_bin_contigs.py <bin_tsv> <fasta_file> <target_bin> <output_fasta>")
        sys.exit(1)

    bin_tsv = sys.argv[1]
    fasta_file = sys.argv[2]
    target_bin = sys.argv[3]
    output_fasta = sys.argv[4]

    extract_contigs(bin_tsv, fasta_file, target_bin, output_fasta)
