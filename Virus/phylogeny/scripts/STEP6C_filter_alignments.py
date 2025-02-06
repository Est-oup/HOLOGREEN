import os
import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import csv
from concurrent.futures import ProcessPoolExecutor

def calculate_gap_percentage(alignment, threshold):
    """
    Calcule le pourcentage de gaps pour chaque position de l'alignement
    et filtre les positions en fonction du seuil donné.
    """
    num_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()

    gap_counts = np.zeros(alignment_length)

    # Calcul du nombre de gaps pour chaque position
    for record in alignment:
        for i, char in enumerate(record.seq):
            if char == '-':
                gap_counts[i] += 1

    # Calcul du pourcentage de gaps
    gap_percentages = gap_counts / num_sequences * 100

    # Filtrage des colonnes
    filtered_alignment = []
    for record in alignment:
        filtered_seq = ''.join([char for i, char in enumerate(record.seq) if gap_percentages[i] < threshold])
        filtered_alignment.append(SeqRecord(Seq(filtered_seq), id=record.id, description=record.description))

    return MultipleSeqAlignment(filtered_alignment), gap_percentages

def process_file(file_info):
    input_file, output_file, threshold = file_info
    alignment = AlignIO.read(input_file, "fasta")
    filtered_alignment, gap_percentages = calculate_gap_percentage(alignment, threshold)

    AlignIO.write(filtered_alignment, output_file, "fasta")

    avg_gap_percentage = np.mean(gap_percentages)
    max_gap_percentage = np.max(gap_percentages)
    min_gap_percentage = np.min(gap_percentages)

    return {
        "file": os.path.basename(input_file),
        "avg_gap_percentage": avg_gap_percentage,
        "max_gap_percentage": max_gap_percentage,
        "min_gap_percentage": min_gap_percentage
    }

def main(input_folder, output_folder, threshold, log_file):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    log_folder = os.path.dirname(log_file)
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)

    file_infos = []
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            input_file = os.path.join(root, file)
            output_file = os.path.join(output_folder, file)
            file_infos.append((input_file, output_file, threshold))

    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        results = executor.map(process_file, file_infos)

    log_data = list(results)

    # Écriture des statistiques dans un fichier CSV
    with open(log_file, 'w', newline='') as csvfile:
        fieldnames = ["file", "avg_gap_percentage", "max_gap_percentage", "min_gap_percentage"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for data in log_data:
            writer.writerow(data)

    print(f"Processing completed. Logs and statistics saved in {log_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python filter_alignment.py <input_folder> <output_folder> <gap_threshold> <log_file>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    gap_threshold = float(sys.argv[3])
    log_file = sys.argv[4]

    main(input_folder, output_folder, gap_threshold, log_file)
