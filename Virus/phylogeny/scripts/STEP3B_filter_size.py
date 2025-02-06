import os
import sys
from Bio import SeqIO

def load_thresholds(threshold_file):
    thresholds = {}
    with open(threshold_file, 'r') as f:
        next(f)  # Ignore header
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                protname, _, threshold = parts
                thresholds[protname] = float(threshold)  # Utiliser float au lieu de int
    return thresholds

def filter_sequences(input_file, output_file, threshold):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if len(record.seq) >= threshold:  # Utiliser >= pour le seuil
                SeqIO.write(record, outfile, "fasta")

def main(input_folder, threshold_file, filtered_output_folder):
    if not os.path.exists(filtered_output_folder):
        os.makedirs(filtered_output_folder)

    thresholds = load_thresholds(threshold_file)

    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith(".fasta"):
                # Gérer le cas où un fichier contient plusieurs fois 'filtered'
                protname = file.replace("_filtered_filtered.fasta", "").replace("_filtered.fasta", "").replace(".fasta", "")
                if protname in thresholds:
                    input_file = os.path.join(root, file)
                    filtered_output_file = os.path.join(filtered_output_folder, file)

                    # Filtrer les séquences
                    filter_sequences(input_file, filtered_output_file, thresholds[protname])
                    print(f"Filtered sequences for {file} with threshold {thresholds[protname]}")
                else:
                    print(f"No threshold found for {file}, skipping.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_and_cluster_orf_size.py <input_folder> <threshold_file> <filtered_output_folder>")
        sys.exit(1)

    input_folder = sys.argv[1]
    threshold_file = sys.argv[2]
    filtered_output_folder = sys.argv[3]

    main(input_folder, threshold_file, filtered_output_folder)
