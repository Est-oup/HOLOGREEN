import os
import sys
from Bio import SeqIO

def process_single_alignment(alignment_path):
    best_hits = {}
    with open(alignment_path, 'r') as f:
        for line in f:
            cols = line.strip().split()
            contig_id = cols[0]
            protein_id = cols[1]
            evalue = float(cols[-2])
            if contig_id not in best_hits or evalue < best_hits[contig_id][1]:
                best_hits[contig_id] = (protein_id, evalue)
    return best_hits

def write_output(best_hits, output_file):
    with open(output_file, 'w') as out_f:
        for contig_id, (protein_id, evalue) in best_hits.items():
            out_f.write(f"{contig_id}\t{protein_id}\t{evalue}\n")



def process_alignments(alignment_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    alignment_files = [os.path.join(alignment_folder, f) for f in os.listdir(alignment_folder) if os.path.isfile(os.path.join(alignment_folder, f)) and not f.startswith('.index')]

    for alignment_file in alignment_files:
        file_prefix = os.path.splitext(os.path.basename(alignment_file))[0]
        best_hits = process_single_alignment(alignment_file)
        output_file = os.path.join(output_folder, f"{file_prefix}.txt")
        write_output(best_hits, output_file)
        print(f"Finished writing output data to {output_file}")



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python STEP1C_process_alignments.py <alignment_folder> <output_folder> <fasta_output_base> <contig_folder> <protein_folder>")
        sys.exit(1)

    alignment_folder = sys.argv[1]
    output_folder = sys.argv[2]


    process_alignments(alignment_folder, output_folder)
