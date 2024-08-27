# -*- coding: utf-8 -*-

import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import threading

lock = threading.Lock()

def run_create_db(input_file, output_folder):
    # Execute la commande mmseqs createdb pour creer une base de donnees.
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    output_file = os.path.join(output_folder, os.path.basename(input_file))
    create_db_command = f'mmseqs createdb {input_file} {output_file}'
    print(f"Executing: {create_db_command}")
    os.system(create_db_command)

def create_blast_db(input_folder, output_folder, max_workers=None):
    # Cree les bases de donnees mmseqs en parallele.
    with lock:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    tasks = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for file_name in os.listdir(input_folder):
            input_file = os.path.join(input_folder, file_name)
            file_output_folder = os.path.join(output_folder, os.path.splitext(file_name)[0])
            tasks.append(executor.submit(run_create_db, input_file, file_output_folder))

        for future in as_completed(tasks):
            future.result()
        print("Base de donnees creee.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python create_db_mmseqs.py <input_folder> <output_folder>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    num_workers = os.cpu_count()
    create_blast_db(input_folder, output_folder, max_workers=num_workers)
