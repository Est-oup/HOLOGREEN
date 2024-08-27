import os
import sys
from concurrent.futures import ProcessPoolExecutor
import threading

lock = threading.Lock()

def run_mmseqs(db_path, ref_path, output_file, output_file_temp, num_threads):
    """
    Exécute la commande mmseqs easy-search pour les fichiers donnés avec le nombre de threads spécifié.
    """
    blast_command = f'mmseqs easy-search {ref_path} {db_path} {output_file} {output_file_temp} --threads {num_threads}'
    print(f"Executing: {blast_command}")
    os.system(blast_command)

def create_output_dir(directory):
    """
    Crée le répertoire s'il n'existe pas déjà.
    """
    with lock:
        if not os.path.exists(directory):
            os.makedirs(directory)

def main(db_query_folder, ref_folder, output_folder, temp_folder):
    """
    Lance les recherches mmseqs en parallèle avec un nombre de threads adapté automatiquement.
    """
    create_output_dir(output_folder)
    create_output_dir(temp_folder)

    # Liste des noms de bases de données et références
    databases = os.listdir(db_query_folder)
    refs = os.listdir(ref_folder)

    # Déterminer le nombre de threads et de processus
    num_threads = os.cpu_count()
    max_workers = num_threads

    # Liste des arguments pour les exécutions parallèles
    tasks = []

    for db_name in databases:
        db_path = os.path.join(db_query_folder, db_name)
        for ref_name in refs:
            ref_path = os.path.join(ref_folder, ref_name)
            output_file = os.path.join(output_folder, f"db_{db_name}_{ref_name}")
            output_file_temp = os.path.join(temp_folder, f"db_{db_name}_{ref_name}")
            tasks.append((db_path, ref_path, output_file, output_file_temp, num_threads))

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
    
    main(db_query_folder, ref_folder, output_folder, temp_folder)
