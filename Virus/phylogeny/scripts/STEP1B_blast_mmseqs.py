import os
import sys
from concurrent.futures import ProcessPoolExecutor
import threading

lock = threading.Lock()

def run_mmseqs(db_path, ref_path, output_file, output_file_temp, num_threads):
    """
    Exécute la commande mmseqs easy-search pour les fichiers donnés avec le nombre de threads spécifié.
    """
    blast_command = f'mmseqs easy-search {db_path} {ref_path} {output_file} {output_file_temp} --threads {num_threads}'
    print(f"Executing command: {blast_command}")
    
    # Exécution de la commande
    result = os.system(blast_command)
    
    # Vérification du résultat
    if result != 0:
        print(f"Error occurred during mmseqs execution: {result}")
    else:
        print(f"mmseqs execution successful for {db_path} vs {ref_path}")

def create_output_dir(directory):
    """
    Crée le répertoire s'il n'existe pas déjà.
    """
    with lock:
        if not os.path.exists(directory):
            print(f"Creating directory: {directory}")
            os.makedirs(directory)
        else:
            print(f"Directory already exists: {directory}")

def extract_marker_name(filename):
    """
    Extrait le nom du marqueur avant l'extension .fasta, avec ou sans '_filtered'.
    """
    print(f"Extracting marker name from: {filename}")
    if '_filtered' in filename:
        return filename.split('_filtered')[0]
    else:
        return filename.split('.fasta')[0]


def main(db_query_folder, ref_folder, output_folder, temp_folder):
    """
    Lance les recherches mmseqs en parallèle avec un nombre de threads adapté automatiquement.
    """
    print(f"Starting process with db_query_folder: {db_query_folder}, ref_folder: {ref_folder}, output_folder: {output_folder}, temp_folder: {temp_folder}")

    # Création des dossiers de sortie
    create_output_dir(output_folder)
    create_output_dir(temp_folder)

    # Liste des noms de bases de données et références
    print("Listing databases and reference files...")
    databases = os.listdir(db_query_folder)
    refs = os.listdir(ref_folder)

    print(f"Found {len(databases)} database files in {db_query_folder}")
    print(f"Found {len(refs)} reference files in {ref_folder}")

    # Déterminer le nombre de threads et de processus
    num_threads = os.cpu_count()
    max_workers = num_threads
    print(f"Using {num_threads} threads")

    # Liste des arguments pour les exécutions parallèles
    tasks = []

    for db_name in databases:
        db_marker = extract_marker_name(db_name)
        db_path = os.path.join(db_query_folder, db_name)
        print(f"Processing database file: {db_name}, marker: {db_marker}")

        for ref_name in refs:
            ref_marker = extract_marker_name(ref_name)
            print(f"Processing reference file: {ref_name}, marker: {ref_marker}")

            # Vérifier si les noms de marqueur correspondent
            if db_marker == ref_marker:
                print(f"Matching markers found: {db_marker}")
                ref_path = os.path.join(ref_folder, ref_name)
                output_file = os.path.join(output_folder, f"{db_name}")
                output_file_temp = os.path.join(temp_folder, f"db_{db_name}_vs_{ref_name}.temp")
                tasks.append((db_path, ref_path, output_file, output_file_temp, num_threads))
            else:
                print(f"No match for {db_marker} and {ref_marker}")

    print(f"Submitting {len(tasks)} tasks for parallel execution")
    
    # Utiliser ProcessPoolExecutor pour paralléliser les tâches
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_mmseqs, *task) for task in tasks]
        for future in futures:
            try:
                future.result()  # Pour lever les exceptions s'il y en a
            except Exception as e:
                print(f"Error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python blast_mmseqs.py <db_query_folder> <ref_folder> <output_folder> <temp_folder>")
        sys.exit(1)

    db_query_folder = sys.argv[1]
    ref_folder = sys.argv[2]
    output_folder = sys.argv[3]
    temp_folder = sys.argv[4]

    print("Script started with the following arguments:")
    print(f"Database query folder: {db_query_folder}")
    print(f"Reference folder: {ref_folder}")
    print(f"Output folder: {output_folder}")
    print(f"Temporary folder: {temp_folder}")
    
    main(db_query_folder, ref_folder, output_folder, temp_folder)
