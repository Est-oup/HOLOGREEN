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

def process_treefile(tree_file, references, output_directory):
    try:
        # Lire et traiter l'arbre IQ-TREE (format treefile)
        tree = Tree(tree_file)
        for leaf in tree:
            # Si le nom n'est pas une séquence de référence, réduire à NOMCONTIG
            if leaf.name not in references:
                parts = leaf.name.split('_')
                if len(parts) > 1:
                    leaf.name = parts[0]
        
        # Écrire le fichier treefile modifié
        output_file = os.path.join(output_directory, os.path.basename(tree_file))
        tree.write(outfile=output_file)
        print(f"Processed {tree_file} -> {output_file}")

    except Exception as e:
        # En cas d'erreur, afficher le fichier problématique et continuer
        print(f"Error processing {tree_file}: {e}")
        
        # Copier le fichier sans modification
        output_file = os.path.join(output_directory, os.path.basename(tree_file))
        with open(tree_file, 'r') as infile, open(output_file, 'w') as outfile:
            outfile.write(infile.read())
        print(f"Copied {tree_file} without modifications to {output_file}")

def main(tree_folder, fasta_directory, output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    # Charger toutes les séquences de référence en mémoire
    references = load_all_reference_sequences(fasta_directory)
    
    for root, _, files in os.walk(tree_folder):
        for file in files:
            # IQ-TREE produit les arbres dans des fichiers .treefile
            if file.endswith(".contree"):
                tree_file = os.path.join(root, file)
                process_treefile(tree_file, references, output_directory)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <tree_folder> <fasta_directory> <output_directory>")
        sys.exit(1)
    
    tree_folder = sys.argv[1]
    fasta_directory = sys.argv[2]
    output_directory = sys.argv[3]

    main(tree_folder, fasta_directory, output_directory)
