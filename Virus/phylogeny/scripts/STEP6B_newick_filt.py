import os
import sys
from Bio import SeqIO
from ete3 import Tree

def load_all_reference_sequences(fasta_directory):
    """
    Lit tous les fichiers FASTA dans un répertoire et retourne un ensemble contenant les noms de toutes les séquences de références.
    """
    references = set()
    for fasta_file in os.listdir(fasta_directory):
        if fasta_file.endswith(".prt"):
            with open(os.path.join(fasta_directory, fasta_file), 'r') as file:
                for line in file:
                    if line.startswith('>'):
                        sequence_name = line.split()[0][1:]  # Enlever le '>'
                        references.add(sequence_name)
    return references

def process_newick_file(newick_file, references, output_directory):
    # Lire et traiter l'arbre Newick
    tree = Tree(newick_file)
    for leaf in tree:
        # Si le nom n'est pas une séquence de référence, réduire à NOMCONTIG
        if leaf.name not in references:
            parts = leaf.name.split('_')
            if len(parts) > 1:
                leaf.name = parts[0]
    
    # Écrire le fichier Newick modifié
    output_file = os.path.join(output_directory, os.path.basename(newick_file))
    tree.write(outfile=output_file)
    print(f"Processed {newick_file} -> {output_file}")

def main(newick_folder, fasta_directory, output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    # Charger toutes les séquences de référence en mémoire
    references = load_all_reference_sequences(fasta_directory)
    
    for root, _, files in os.walk(newick_folder):
        for file in files:
            if file.endswith(".nw") or file.endswith(".newick"):
                newick_file = os.path.join(root, file)
                process_newick_file(newick_file, references, output_directory)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <newick_folder> <fasta_directory> <output_directory>")
        sys.exit(1)
    
    newick_folder = sys.argv[1]
    fasta_directory = sys.argv[2]
    output_directory = sys.argv[3]

    main(newick_folder, fasta_directory, output_directory)
