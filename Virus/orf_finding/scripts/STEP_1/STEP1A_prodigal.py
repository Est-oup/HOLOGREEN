import os
import sys
import subprocess
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO

def split_fasta(input_file, output_dir, chunk_size=1000):
    """
    Divise un fichier FASTA en plusieurs segments.
    """
    os.makedirs(output_dir, exist_ok=True)
    chunk_files = []
    try:
        with open(input_file, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            for i in range(0, len(records), chunk_size):
                chunk_records = records[i:i + chunk_size]
                chunk_file = os.path.join(output_dir, f"chunk_{i//chunk_size}.fasta")
                with open(chunk_file, "w") as chunk_handle:
                    SeqIO.write(chunk_records, chunk_handle, "fasta")
                chunk_files.append(chunk_file)
    except Exception as e:
        print(f"Error while splitting fasta: {e}")
    return chunk_files

def run_prodigal(input_file, output_file):
    """
    Exécute Prodigal sur un fichier d'entrée et écrit le résultat dans un fichier de sortie.
    """
    command = f'prodigal -i {input_file} -a {output_file} -p meta'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        return f"Prédiction d'ORF créée pour le fichier {os.path.basename(input_file)} dans {output_file}"
    else:
        return f"Erreur lors de la prédiction d'ORF pour le fichier {os.path.basename(input_file)}: {result.stderr}"

def count_orfs(fasta_file):
    """
    Compte le nombre d'ORF dans un fichier FASTA.
    """
    count = 0
    try:
        with open(fasta_file, "r") as handle:
            for _ in SeqIO.parse(handle, "fasta"):
                count += 1
    except Exception as e:
        print(f"Error while counting ORFs: {e}")
    return count

def prodigalsearch(input_folder, output_folder, temp_folder, log_file, chunk_size=1000, max_workers=None):
    """
    Parcours de chaque fichier dans le dossier d'entrée, divise les fichiers, exécute Prodigal en parallèle,
    et fusionne les résultats.
    """
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(temp_folder, exist_ok=True)

    input_files = [os.path.join(input_folder, file_name) for file_name in os.listdir(input_folder) if os.path.isfile(os.path.join(input_folder, file_name))]

    if not input_files:
        print("Aucun fichier d'entrée trouvé.")
        return

    total_orfs = 0

    for input_file in input_files:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        chunk_files = split_fasta(input_file, temp_folder, chunk_size)
        output_files = [os.path.join(temp_folder, f"{base_name}_chunk_{i}.fasta") for i in range(len(chunk_files))]
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(run_prodigal, chunk, output) for chunk, output in zip(chunk_files, output_files)]
            for future in futures:
                print(future.result())
        
        # Fusion des fichiers de sortie en un seul fichier de sortie final
        final_output_file = os.path.join(output_folder, f"{base_name}.fasta")
        with open(final_output_file, "w") as outfile:
            for output_file in output_files:
                with open(output_file, "r") as infile:
                    outfile.write(infile.read())
        
        # Compter les ORF dans le fichier de sortie final
        orf_count = count_orfs(final_output_file)
        total_orfs += orf_count
        print(f"Prédiction d'ORF fusionnée pour le fichier {input_file} dans {final_output_file} avec {orf_count} ORFs trouvés.")

    # Écrire le nombre total d'ORF trouvés dans le fichier log
    with open(log_file, "w") as log_handle:
        log_handle.write(f"Nombre total d'ORFs trouvés: {total_orfs}\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <input_folder> <output_folder> <temp_folder> <log_file>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    temp_folder = sys.argv[3]  # dossier temporaire pour stocker les segments
    log_file = sys.argv[4]     # fichier log pour enregistrer le nombre d'ORF trouvés

    # Utiliser tous les CPU disponibles
    num_workers = os.cpu_count()

    prodigalsearch(input_folder, output_folder, temp_folder, log_file, max_workers=num_workers)

