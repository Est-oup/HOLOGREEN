import os
import sys
import pandas as pd
from Bio import SeqIO

def calculate_min_sequence_lengths(input_folder, output_file):
    data = []

    for root, dirs, files in os.walk(input_folder):
        for file in files:
            file_path = os.path.join(root, file)
            min_length = float('inf')
            for record in SeqIO.parse(file_path, "fasta"):
                seq_length = len(record.seq)
                if seq_length < min_length:
                    min_length = seq_length
                
            protein_name = os.path.splitext(file)[0]
            data.append({
                "Protein": protein_name,
                "Min_Length": min_length,
                "0.75_Min_Length": 0.75 * min_length
            })
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Save data to TSV (Tab-separated values)
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False, sep='\t')

    print("Analysis completed. Data saved in {}".format(output_file))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python calculate_min_lengths.py <input_folder> <output_file>")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    output_file = sys.argv[2]

    calculate_min_sequence_lengths(input_folder, output_file)
