import os
import sys
import subprocess

def create_output_dirs(output_aln_dir, temp_aln_dir):
    os.makedirs(output_aln_dir, exist_ok=True)
    os.makedirs(temp_aln_dir, exist_ok=True)

def run_mmseqs(db_path, fasta_path, output_file, temp_file, num_threads=16):
    """
    Exécute la commande mmseqs easy-search pour les fichiers donnés avec le nombre de threads spécifié.
    """
    command = [
        'mmseqs', 'easy-search',
        fasta_path, db_path,
        output_file, temp_file,
        '--threads', str(num_threads)
    ]
    print(f"Executing: {' '.join(command)}")
    subprocess.run(command, check=True)

def find_db_file(db_dir, db_name):
    # Vérifier les extensions possibles
    possible_extensions = ["", ".dbtype", ".index", ".lookup", ".source", "_h", "_h.dbtype", "_h.index"]
    for ext in possible_extensions:
        db_path = os.path.join(db_dir, db_name + ext)
        if os.path.exists(db_path):
            return db_path
    return None

def main(db_dir, fasta_dir, output_aln_dir, temp_aln_dir):
    # Create output directories
    create_output_dirs(output_aln_dir, temp_aln_dir)
    
    for fasta_file in os.listdir(fasta_dir):
        if fasta_file.endswith(".fasta"):
            fasta_path = os.path.join(fasta_dir, fasta_file)
            base_name = fasta_file.split(".fasta")[0]

            # Extract base name to find corresponding database
            db_key = "_".join(base_name.split("_")[:-1])
            db_name = db_key.replace("_", "-")

            # Trouver le fichier de base de données
            db_path = find_db_file(os.path.join(db_dir, db_name), db_name)

            # Debugging output
            print(f"Processing {fasta_file}")
            print(f"Constructed db_key: {db_key}")
            print(f"Expected db_path: {db_path}")

            if db_path:
                output_file = os.path.join(output_aln_dir, f"{base_name}_filtered")
                temp_file = os.path.join(temp_aln_dir, f"{base_name}_filtered.tmp")

                run_mmseqs(db_path, fasta_path, output_file, temp_file)
            else:
                print(f"Database path for {fasta_file} does not exist")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <db_dir> <fasta_dir> <output_aln_dir> <temp_aln_dir>")
        sys.exit(1)
    
    db_dir = sys.argv[1]
    fasta_dir = sys.argv[2]
    output_aln_dir = sys.argv[3]
    temp_aln_dir = sys.argv[4]
    
    main(db_dir, fasta_dir, output_aln_dir, temp_aln_dir)
