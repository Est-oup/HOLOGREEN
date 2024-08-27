import os
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor

def align_sequences(fasta_file, output_folder):
    file_basename = os.path.basename(fasta_file).replace(".fasta", "")
    aligned_file = os.path.join(output_folder, f"{file_basename}.fasta")
    
    mafft_cmd = f"mafft --anysymbol {fasta_file} > {aligned_file}"
    
    try:
        subprocess.run(mafft_cmd, shell=True, check=True)
        if os.path.exists(aligned_file):
            print(f"Alignment completed for {fasta_file}")
        else:
            print(f"Error: Aligned file {aligned_file} not found after MAFFT execution.")
            return
    except subprocess.CalledProcessError as e:
        print(f"Error during alignment with MAFFT for {fasta_file}: {e}", file=sys.stderr)
        return

def main(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    fasta_files = []
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith(".fasta"):
                fasta_files.append(os.path.join(root, file))

    print(f"Found {len(fasta_files)} FASTA files to process.")
    
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        executor.map(lambda f: align_sequences(f, output_folder), fasta_files)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python STEP1E_align_sequences.py <input_folder> <output_folder>")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    print(f"Starting alignment process with the following parameters:\nInput folder: {input_folder}\nOutput folder: {output_folder}")
    main(input_folder, output_folder)
