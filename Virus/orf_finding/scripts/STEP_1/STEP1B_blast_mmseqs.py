import os
import sys
from concurrent.futures import ProcessPoolExecutor
import threading

lock = threading.Lock()

def run_mmseqs(query_file, db_path, output_file, output_file_temp, num_threads):
    """
    Execute the mmseqs easy-search command for the given files with the specified number of threads.
    """
    blast_command = f'mmseqs easy-search {query_file} {db_path} {output_file} {output_file_temp} --threads {num_threads}'
    print(f"Executing: {blast_command}")
    os.system(blast_command)

def create_output_dir(directory):
    """
    Create the directory if it doesn't already exist.
    """
    with lock:
        if not os.path.exists(directory):
            os.makedirs(directory)

def main(query_path, db_folder, output_folder, temp_folder):
    """
    Launch mmseqs searches in parallel with an automatically adjusted number of threads.
    """
    create_output_dir(output_folder)
    create_output_dir(temp_folder)

    # List database paths
    db_paths = [os.path.join(db_folder, db_name, db_name) for db_name in os.listdir(db_folder) if os.path.isdir(os.path.join(db_folder, db_name))]

    # Determine the number of threads and processes
    num_threads = os.cpu_count() // len(db_paths)
    if num_threads < 1:
        num_threads = 1

    # List of arguments for parallel executions
    tasks = []

    for db_path in db_paths:
        db_name = os.path.basename(os.path.dirname(db_path))  # Extract the DB name
        output_file = os.path.join(output_folder, f"db_{db_name}_{os.path.basename(query_path)}")
        output_file_temp = os.path.join(temp_folder, f"db_{db_name}_{os.path.basename(query_path)}")
        tasks.append((query_path, db_path, output_file, output_file_temp, num_threads))

    # Use ProcessPoolExecutor to parallelize tasks
    with ProcessPoolExecutor(max_workers=len(tasks)) as executor:
        futures = [executor.submit(run_mmseqs, *task) for task in tasks]
        for future in futures:
            try:
                future.result()  # To raise exceptions if any
            except Exception as e:
                print(f"Error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python STEP1B_blast_mmseqs.py <query_path> <db_folder> <output_folder> <temp_folder>")
        sys.exit(1)

    query_path = sys.argv[1]
    db_folder = sys.argv[2]
    output_folder = sys.argv[3]
    temp_folder = sys.argv[4]
    
    main(query_path, db_folder, output_folder, temp_folder)
