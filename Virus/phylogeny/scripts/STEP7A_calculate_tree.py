import os
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor

def build_tree(aligned_file, output_folder):
    file_basename = os.path.basename(aligned_file).replace(".fasta", "")
    tree_output_dir = os.path.join(output_folder, file_basename)
    
    # Création du sous-dossier pour les fichiers de sortie IQ-TREE
    if not os.path.exists(tree_output_dir):
        os.makedirs(tree_output_dir)
    
    # Commande IQ-TREE avec les options et le préfixe de sortie dans le sous-dossier
    iqtree_cmd = f"iqtree2 -s {aligned_file} -m MFP -msub viral -B 1000 -nt AUTO -pre {os.path.join(tree_output_dir, file_basename)}"

    try:
        # Exécution de la commande IQ-TREE
        subprocess.run(iqtree_cmd, shell=True, check=True)
        
        print(f"Tree and associated files generated for {aligned_file} in {tree_output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"Error during tree building with IQ-TREE for {aligned_file}: {e}", file=sys.stderr)

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
        print("Usage: python build_phylo_trees_with_iqtree.py <input_folder> <output_folder>")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    print(f"Starting tree building process with the following parameters:\nInput folder: {input_folder}\nOutput folder: {output_folder}")
    main(input_folder, output_folder)
