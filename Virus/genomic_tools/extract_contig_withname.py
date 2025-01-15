import sys
from Bio import SeqIO

def extract_sequences_by_pattern(input_fasta, pattern, output_fasta):
    """
    Extrait les séquences d'un fichier FASTA dont l'en-tête commence par un motif donné.

    :param input_fasta: Chemin vers le fichier FASTA d'entrée
    :param pattern: Motif à rechercher en début d'en-tête
    :param output_fasta: Chemin vers le fichier FASTA de sortie
    """
    try:
        # Lecture des séquences du fichier FASTA d'entrée
        with open(input_fasta, "r") as infile:
            sequences = list(SeqIO.parse(infile, "fasta"))
        
        # Filtrage des séquences dont l'en-tête commence par le motif donné
        matching_sequences = [seq for seq in sequences if seq.id.startswith(pattern)]
        
        # Écriture des séquences correspondantes dans le fichier FASTA de sortie
        with open(output_fasta, "w") as outfile:
            SeqIO.write(matching_sequences, outfile, "fasta")
        
        print(f"Extraction terminée : {len(matching_sequences)} séquences extraites et écrites dans {output_fasta}.")
    except Exception as e:
        print(f"Erreur lors du traitement : {e}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage : python extract_fasta_by_pattern.py <input_fasta> <pattern> <output_fasta>")
    else:
        input_fasta = sys.argv[1]
        pattern = sys.argv[2]
        output_fasta = sys.argv[3]
        extract_sequences_by_pattern(input_fasta, pattern, output_fasta)
