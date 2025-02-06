import os
import sys
from Bio import SeqIO
from collections import defaultdict

def parse_gff(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("##STA"):
                sequences.append(line.strip().split('\t')[1])
    return sequences

def write_fasta(sequences, output_file, protein_type):
    with open(output_file, 'w') as file:
        for i, (seq, header) in enumerate(sequences):
            file.write(f">{header}_{protein_type}_{i+1}\n")  # Inclut le nom de la protéine dans l'identifiant
            file.write(f"{seq}\n")

def write_combined_fasta(sequences, output_file):
    with open(output_file, 'w') as file:
        for i, (seq, header) in enumerate(sequences):
            contig_name = header.split("_")[0]  # Extrait le nom du contig
            file.write(f">{contig_name}_{i+1}\n")
            file.write(f"{seq}\n")

def process_files(input_dir, output_dir, combined_file, log_file):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    total_files_processed = 0
    total_orfs_found = 0
    combined_sequences = []

    protein_sequences = defaultdict(list)

    with open(log_file, 'w') as log:
        log.write("File\tORFs Found\n")

        for root, dirs, files in os.walk(input_dir):
            for name in files:
                if name.endswith(".gff"):
                    file_path = os.path.join(root, name)
                    sequences = parse_gff(file_path)
                    if sequences:
                        protein_type = os.path.basename(root)  # Récupère le nom du dossier comme type de protéine
                        protein_output_dir = os.path.join(output_dir, protein_type)
                        if not os.path.exists(protein_output_dir):
                            os.makedirs(protein_output_dir)

                        output_file = os.path.join(protein_output_dir, f"{os.path.splitext(name)[0]}.fasta")
                        headers = [f"{os.path.splitext(name)[0]}_{protein_type}_{i+1}" for i in range(len(sequences))]
                        write_fasta(list(zip(sequences, headers)), output_file, protein_type)

                        protein_sequences[protein_type].extend(zip(sequences, headers))
                        combined_sequences.extend(zip(sequences, headers))

                        log.write(f"{file_path}\t{len(sequences)}\n")
                        print(f"Processed {file_path}, found {len(sequences)} ORFs.")

                        total_files_processed += 1
                        total_orfs_found += len(sequences)

        for protein_type, seqs in protein_sequences.items():
            protein_output_file = os.path.join(output_dir, f"{protein_type}.fasta")
            write_fasta(seqs, protein_output_file, protein_type)

        write_combined_fasta(combined_sequences, combined_file)

        log.write(f"\nSummary:\n")
        log.write(f"Total files processed: {total_files_processed}\n")
        log.write(f"Total ORFs found: {total_orfs_found}\n")

def main():
    if len(sys.argv) != 5:
        print("Usage: python STEP1E_extract_miniprot_results.py <input_dir> <output_dir> <combined_file> <log_file>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    combined_file = sys.argv[3]
    log_file = sys.argv[4]

    process_files(input_dir, output_dir, combined_file, log_file)

if __name__ == "__main__":
    main()
