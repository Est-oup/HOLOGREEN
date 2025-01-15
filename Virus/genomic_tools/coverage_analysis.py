import sys
import pandas as pd

def parse_sam(sam_file):
    """
    Parse a SAM file and calculate alignment coverage for each reference.
    """
    coverage = {}
    with open(sam_file, 'r') as file:
        for line in file:
            if line.startswith("@"):  # Skip header lines
                continue
            fields = line.split("\t")
            ref_name = fields[2]  # Reference name
            ref_start = int(fields[3])  # Start position (1-based)
            cigar = fields[5]  # CIGAR string

            # Calculate alignment length from CIGAR
            alignment_length = sum([int(num) for num in "".join([c if c.isdigit() else " " for c in cigar]).split()])

            if ref_name not in coverage:
                coverage[ref_name] = 0
            coverage[ref_name] += alignment_length

    return coverage

def calculate_coverage(reference_fasta, coverage_dict):
    """
    Compare alignment coverage to the total length of untrimmed contigs.
    """
    ref_lengths = {}
    with open(reference_fasta, 'r') as file:
        seq_id = None
        seq_length = 0
        for line in file:
            if line.startswith(">"):
                if seq_id:
                    ref_lengths[seq_id] = seq_length
                seq_id = line[1:].strip()
                seq_length = 0
            else:
                seq_length += len(line.strip())
        if seq_id:
            ref_lengths[seq_id] = seq_length

    results = []
    for seq_id, total_length in ref_lengths.items():
        aligned_length = coverage_dict.get(seq_id, 0)
        coverage_percentage = (aligned_length / total_length) * 100
        results.append({
            "Sequence_ID": seq_id,
            "Total_Length": total_length,
            "Aligned_Length": aligned_length,
            "Coverage_Percentage": coverage_percentage
        })

    return pd.DataFrame(results)

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <reference_fasta> <trimmed_fasta>")
        sys.exit(1)

    reference_fasta = sys.argv[1]  # Untrimmed contigs
    trimmed_fasta = sys.argv[2]  # Trimmed contigs

    # Step 1: Run Minimap2 alignment (you need minimap2 installed in your system)
    sam_file = "alignment.sam"
    command = f"minimap2 -x asm5 -a {reference_fasta} {trimmed_fasta} > {sam_file}"
    print(f"Running: {command}")
    import os
    os.system(command)

    # Step 2: Parse SAM file
    coverage_dict = parse_sam(sam_file)

    # Step 3: Calculate coverage
    coverage_results = calculate_coverage(reference_fasta, coverage_dict)

    # Step 4: Display results
    print("Coverage Analysis Results:")
    print(coverage_results)

    # Optionally, save to a file
    output_file = "coverage_results.csv"
    coverage_results.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
