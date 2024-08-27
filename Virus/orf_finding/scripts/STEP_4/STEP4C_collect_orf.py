# -*- coding: utf-8 -*-

import os
import sys
from Bio import SeqIO
from collections import defaultdict

def list_files(directory):
    # Liste tous les fichiers dans le repertoire donne
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

def parse_protein_files(files):
    # Regroupe les fichiers par nom de proteine, sans l'annotation specifique
    protein_files = defaultdict(list)
    for file in files:
        # Extrait le nom de la proteine en retirant l'annotation et l'extension
        base_name = file.replace('db_ORF_pred_derep_', '').replace('_id_sequences.fasta', '')
        protein_name = base_name.split('-')[0]
        protein_files[protein_name].append(file)
    return protein_files

def merge_fasta_files(input_dir, output_dir, protein_name, files):
    # Fusionne les fichiers fasta pour une proteine donnee, sans doublons
    sequences = {}
    for file in files:
        file_path = os.path.join(input_dir, file)
        for record in SeqIO.parse(file_path, "fasta"):
            sequences[str(record.seq)] = record.description
    
    output_file = os.path.join(output_dir, f"{protein_name}.fasta")
    with open(output_file, 'w') as output_handle:
        for seq, description in sequences.items():
            output_handle.write(f">{description}\n{seq}\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_directory> <output_directory>")
        sys.exit(1)
    
    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    files = list_files(input_directory)
    protein_files = parse_protein_files(files)
    
    for protein_name, file_list in protein_files.items():
        if len(file_list) > 1:
            merge_fasta_files(input_directory, output_directory, protein_name, file_list)

if __name__ == "__main__":
    main()

