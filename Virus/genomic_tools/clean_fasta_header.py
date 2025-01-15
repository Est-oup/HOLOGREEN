import sys
import os

def clean_fasta_headers(input_fasta, output_fasta):
    """Nettoie les en-têtes d'un fichier FASTA en gardant uniquement la première partie de l'en-tête."""
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith(">"):
                header = line.split()[0]
                outfile.write(header + "\n")
            else:
                outfile.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python clean_fasta.py <input_fasta> <output_fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]

    if not os.path.exists(input_fasta):
        print(f"Le fichier {input_fasta} n'existe pas.")
        sys.exit(1)

    clean_fasta_headers(input_fasta, output_fasta)
    print(f"Les en-têtes du fichier FASTA ont été nettoyées et enregistrées dans {output_fasta}.")
