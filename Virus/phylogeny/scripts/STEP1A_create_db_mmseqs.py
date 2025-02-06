import os
import sys
import subprocess

def run_create_db(input_file, output_file):
    """
    Exécute la commande mmseqs createdb pour chaque fichier séquence.
    """
    create_db_command = f'mmseqs createdb {input_file} {output_file}'
    print(f"Executing: {create_db_command}")
    
    # Capture la sortie et les erreurs
    process = subprocess.Popen(create_db_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    # Afficher les sorties
    if stdout:
        print(f"Output: {stdout.decode()}")
    if stderr:
        print(f"Error: {stderr.decode()}")
    if process.returncode != 0:
        print(f"Command failed with return code {process.returncode}")

def create_blast_db(input_folder, output_folder):
    """
    Parcourt chaque fichier dans le dossier d'entrée et crée une base de données pour chacun.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for file_name in os.listdir(input_folder):
        input_file = os.path.join(input_folder, file_name)
        if os.path.isfile(input_file):  # Vérifie que c'est bien un fichier
            output_file = os.path.join(output_folder, os.path.splitext(file_name)[0])
            run_create_db(input_file, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python create_db_mmseqs.py <input_folder> <output_folder>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    create_blast_db(input_folder, output_folder)
