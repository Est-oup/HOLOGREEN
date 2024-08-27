import os
import sys
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed

def rename_sequences(file_path, new_header):
    renamed_sequences = []
    with open(file_path, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            record.description = new_header
            record.id = new_header.split(" ")[0]
            renamed_sequences.append(record)
    return renamed_sequences

def process_folder_or_file(path, new_header):
    sequences = []
    if os.path.isdir(path):
        with ThreadPoolExecutor() as executor:
            futures = [executor.submit(rename_sequences, os.path.join(path, file_name), new_header) 
                       for file_name in os.listdir(path)]
            for future in as_completed(futures):
                sequences.extend(future.result())
    else:
        sequences = rename_sequences(path, new_header)
    return sequences

def merge_sequences(sequences1, sequences2, output_file):
    merged_sequences = sequences1 + sequences2
    with open(output_file, "w") as output_handle:
        SeqIO.write(merged_sequences, output_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py bact_DNA_pol_folder ref_folder output_folder")
        sys.exit(1)

    bact_DNA_pol_folder = sys.argv[1]
    ref_folder = sys.argv[2]
    output_folder = sys.argv[3]

    # Renommer les séquences du dossier bact_DNA_pol_folder avec le format ">Bacteria"
    bact_sequences = process_folder_or_file(bact_DNA_pol_folder, "Bacteria")

    # Renommer les séquences du dossier ref_folder avec le format ">Virus"
    ref_sequences = process_folder_or_file(ref_folder, "Virus")

    # Fusionner les séquences renommées en un seul fichier
    output_file_path = os.path.join(output_folder, "merged_sequences.fasta")
    merge_sequences(bact_sequences, ref_sequences, output_file_path)

    print(f"Merged sequences written to {output_file_path}")
