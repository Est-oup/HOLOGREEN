import os
import sys
import subprocess
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

def cluster_sequences(input_file, output_fasta, output_clstr_folder, identity_threshold, log_file):
    cdhit_cmd = f"cd-hit -i {input_file} -o {output_fasta} -c {identity_threshold} -d 0"
    with open(log_file, 'w') as logf:
        try:
            subprocess.run(cdhit_cmd, shell=True, check=True, stdout=logf, stderr=subprocess.STDOUT)
            output_clstr = f"{output_fasta}.clstr"
            final_clstr_path = os.path.join(output_clstr_folder, os.path.basename(output_clstr))
            if os.path.exists(output_fasta) and os.path.exists(output_clstr):
                os.rename(output_clstr, final_clstr_path)
                print(f"Clustering completed for {input_file}")
            else:
                print(f"Error: Clustering output files not found for {input_file}")
                with open(log_file, 'r') as lf:
                    print(lf.read())
        except subprocess.CalledProcessError as e:
            print(f"Error during clustering with CD-HIT for {input_file}: {e}", file=sys.stderr)
            with open(log_file, 'r') as lf:
                print(lf.read())

def main(input_folder, threshold_file, filtered_output_folder, clustered_output_folder, clstr_output_folder, log_folder, identity_threshold):
    if not os.path.exists(filtered_output_folder):
        os.makedirs(filtered_output_folder)
    if not os.path.exists(clustered_output_folder):
        os.makedirs(clustered_output_folder)
    if not os.path.exists(clstr_output_folder):
        os.makedirs(clstr_output_folder)
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)

    thresholds = load_thresholds(threshold_file)

    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith(".fasta"):
                protname = file.replace("_filtered.fasta", "").replace(".fasta", "")
                if protname in thresholds:
                    input_file = os.path.join(root, file)
                    filtered_output_file = os.path.join(filtered_output_folder, file)
                    clustered_output_file = os.path.join(clustered_output_folder, file)
                    log_file = os.path.join(log_folder, f"{file}.log")

                    # Filtrer les séquences
                    filter_sequences(input_file, filtered_output_file, thresholds[protname])
                    print(f"Filtered sequences for {file} with threshold {thresholds[protname]}")

                    # Clusteriser les séquences
                    cluster_sequences(filtered_output_file, clustered_output_file, clstr_output_folder, identity_threshold, log_file)
                else:
                    print(f"No threshold found for {file}, skipping.")

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python filter_and_cluster_orf_size.py <input_folder> <threshold_file> <filtered_output_folder> <clustered_output_folder> <clstr_output_folder> <log_folder> <identity_threshold>")
        sys.exit(1)

    input_folder = sys.argv[1]
    threshold_file = sys.argv[2]
    filtered_output_folder = sys.argv[3]
    clustered_output_folder = sys.argv[4]
    clstr_output_folder = sys.argv[5]
    log_folder = sys.argv[6]
    identity_threshold = float(sys.argv[7])

    main(input_folder, threshold_file, filtered_output_folder, clustered_output_folder, clstr_output_folder, log_folder, identity_threshold)
