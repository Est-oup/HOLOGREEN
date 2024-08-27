import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import threading

lock = threading.Lock()

def run_create_db(input_file, output_file):
    """
    Exécute la commande mmseqs createdb pour créer une base de données.
    """
    create_db_command = f'mmseqs createdb {input_file} {output_file}'
    print(f"Executing: {create_db_command}")
    os.system(create_db_command)

def create_blast_db(input_file, output_folder, max_workers=None):
    """
    Crée la base de données mmseqs pour le fichier d'entrée.
    """
    with lock:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
    
    tasks = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        output_file = os.path.join(output_folder, os.path.splitext(os.path.basename(input_file))[0])
        tasks.append(executor.submit(run_create_db, input_file, output_file))
        
        for future in as_completed(tasks):
            future.result()
            print("Base de données créée.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python create_db_mmseqs.py <input_file> <output_folder>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_folder = sys.argv[2]
    
    num_workers = os.cpu_count()
    create_blast_db(input_file, output_folder, max_workers=num_workers)
