import os
import sys
from Bio import SeqIO

def replace_stars_in_fasta(input_fasta, output_fasta, replace_char):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if not line.startswith('>'):
                line = line.replace('*', replace_char)
            outfile.write(line)

def compile_fasta_files(folder1, folder2, output_file):
    all_sequences = []
    total_sequences_folder1 = 0
    total_sequences_folder2 = 0

    # Parcourir le premier dossier
    for file_name in os.listdir(folder1):
        if file_name.endswith(".fasta"):
            file_path = os.path.join(folder1, file_name)
            sequences = list(SeqIO.parse(file_path, "fasta"))
            total_sequences_folder1 += len(sequences)
            all_sequences.extend(sequences)

    # Parcourir le deuxième dossier
    for file_name in os.listdir(folder2):
        if file_name.endswith(".fasta"):
            file_path = os.path.join(folder2, file_name)
            sequences = list(SeqIO.parse(file_path, "fasta"))
            total_sequences_folder2 += len(sequences)
            all_sequences.extend(sequences)

    # Écrire toutes les séquences dans le fichier de sortie
    with open(output_file, "w") as output_handle:
        SeqIO.write(all_sequences, output_handle, "fasta")

    return total_sequences_folder1, total_sequences_folder2, len(all_sequences)

def dereplicate(input_file, output_file, replace_char):
    log_dir = os.path.join(os.path.dirname(output_file), "logs")
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    log_file = os.path.join(log_dir, os.path.basename(output_file) + ".log")
    clstr_file = os.path.join(log_dir, os.path.basename(output_file) + ".clstr")

    # Remplacer les étoiles par "X" dans le fichier d'entrée FASTA
    modified_input_file = input_file + ".modified"
    replace_stars_in_fasta(input_file, modified_input_file, replace_char)

    # Exécuter CD-HIT sur le fichier modifié
    os.system(f"cd-hit-est -i {modified_input_file} -o {output_file} -c 0.99 -M 0 -d 0 > {log_file}")
    os.rename(output_file + ".clstr", clstr_file)

    # Remplacer les "X" par des étoiles dans le fichier de sortie FASTA
    final_output_file = output_file + ".final"
    replace_stars_in_fasta(output_file, final_output_file, replace_char)

    # Compter les séquences dédupliquées
    dereplicated_sequences = list(SeqIO.parse(final_output_file, "fasta"))
    return len(dereplicated_sequences)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script.py folder1 folder2 output_file dereplicated_output_file super_log_file")
        sys.exit(1)

    folder1 = sys.argv[1]
    folder2 = sys.argv[2]
    output_file = sys.argv[3]
    dereplicated_output_file = sys.argv[4]
    super_log_file = sys.argv[5]

    replace_char = "X"

    total_seq1, total_seq2, total_combined = compile_fasta_files(folder1, folder2, output_file)
    total_dereplicated = dereplicate(output_file, dereplicated_output_file, replace_char)

    with open(super_log_file, "w") as log_handle:
        log_handle.write(f"Total sequences from {folder1}: {total_seq1}\n")
        log_handle.write(f"Total sequences from {folder2}: {total_seq2}\n")
        log_handle.write(f"Total combined sequences: {total_combined}\n")
        log_handle.write(f"Total dereplicated sequences: {total_dereplicated}\n")
        log_handle.write(f"Reduction: {total_combined - total_dereplicated} sequences removed\n")

