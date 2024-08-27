import os
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor

def build_tree(aligned_file, output_folder):
    file_basename = os.path.basename(aligned_file).replace(".fasta", "")
    tree_file = os.path.join(output_folder, f"{file_basename}_tree.newick")
    
    fasttree_cmd = f"FastTree {aligned_file} > {tree_file}"
    
    try:
        subprocess.run(fasttree_cmd, shell=True, check=True)
        if os.path.exists(tree_file):
            print(f"Tree built for {aligned_file}")
        else:
            print(f"Error: Tree file {tree_file} not found after FastTree execution.")
    except subprocess.CalledProcessError as e:
        print(f"Error during tree building with FastTree for {aligned_file}: {e}", file=sys.stderr)

def main(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    aligned_files = []
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith(".fasta"):
                aligned_files.append(os.path.join(root, file))

    print(f"Found {len(aligned_files)} aligned FASTA files to process.")
    
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        executor.map(lambda f: build_tree(f, output_folder), aligned_files)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python build_phylo_trees_with_fasttree.py <input_folder> <output_folder>")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    print(f"Starting tree building process with the following parameters:\nInput folder: {input_folder}\nOutput folder: {output_folder}")
    main(input_folder, output_folder)
