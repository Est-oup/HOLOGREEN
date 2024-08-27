import os
import sys
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

def filter_and_copy_file(src_path, dst_path):
    """
    Copie un fichier de 'src_path' vers 'dst_path'.
    """
    shutil.copy(src_path, dst_path)
    return f"Copied {src_path} to {dst_path}"

def filter_and_copy_files(input_folder, output_folder, max_workers=8):
    """
    Copie les fichiers de 'input_folder' qui contiennent 'dna' et 'pol' dans le nom,
    indépendamment de la casse, vers 'output_folder' en utilisant la parallélisation.
    """
    os.makedirs(output_folder, exist_ok=True)

    tasks = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for filename in os.listdir(input_folder):
            if 'dna' in filename.lower() and 'pol' in filename.lower():
                src_path = os.path.join(input_folder, filename)
                dst_path = os.path.join(output_folder, filename)
                tasks.append(executor.submit(filter_and_copy_file, src_path, dst_path))

        for future in as_completed(tasks):
            print(future.result())

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_folder> <output_folder>")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    # Utiliser le nombre de threads égal au nombre de CPU pour maximiser la performance
    num_workers = os.cpu_count()

    filter_and_copy_files(input_folder, output_folder, max_workers=num_workers)
