import os
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

def create_directory_if_not_exists(directory):
    """Crée un répertoire s'il n'existe pas déjà."""
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Répertoire créé : {directory}")
    else:
        print(f"Répertoire déjà existant : {directory}")

def create_blast_db(input_file_path, db_path):
    """Crée une base de données MMseqs2."""
    create_directory_if_not_exists(os.path.dirname(db_path))

    create_db_cmd = f"mmseqs createdb {input_file_path} {db_path}"
    try:
        subprocess.run(create_db_cmd, shell=True, check=True)
        print(f"Base de données créée pour {input_file_path} à {db_path}")
    except subprocess.CalledProcessError as e:
        print(f"Erreur lors de la création de la base de données pour {input_file_path}: {e}")

def run_mmseqs_search(input_file_path, db_path, aln_file, tmp_file, max_workers):
    """Exécute une recherche MMseqs2."""
    create_directory_if_not_exists(os.path.dirname(aln_file))
    create_directory_if_not_exists(os.path.dirname(tmp_file))

    search_cmd = (
        f"mmseqs easy-search {input_file_path} {db_path} {aln_file} {tmp_file} "
        f"--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits --threads {max_workers}"
    )
    try:
        subprocess.run(search_cmd, shell=True, check=True)
        print(f"Recherche terminée pour {input_file_path}")
    except subprocess.CalledProcessError as e:
        print(f"Erreur lors de la recherche MMseqs2 pour {input_file_path}: {e}")

def process_file(filename, input_files_folder, db_folder, aln_folder, result_folder, tmp_folder, max_workers):
    """Processus complet pour un fichier donné."""
    base_name = os.path.splitext(filename)[0]
    input_file_path = os.path.join(input_files_folder, filename)
    db_path = os.path.join(db_folder, base_name)
    aln_file = os.path.join(aln_folder, base_name)
    tmp_file = os.path.join(tmp_folder, base_name)

    # Créer la base de données
    create_blast_db(input_file_path, db_path)

    # Exécuter la recherche MMseqs pour chaque fichier
    run_mmseqs_search(input_file_path, db_path, aln_file, tmp_file, max_workers)

    # Analyser les résultats et enregistrer les pires scores
    worst_scores = extract_worst_scores(aln_file)
    save_worst_scores(worst_scores, os.path.join(result_folder, f"{base_name}_worst_scores.txt"))

def extract_worst_scores(aln_file):
    """Extrait les pires scores d'alignement."""
    worst_scores = {}
    try:
        with open(aln_file, 'r') as f:
            for line in f:
                parts = line.split()
                seq_id = parts[0]
                score = float(parts[11])

                if seq_id not in worst_scores or score < worst_scores[seq_id]:
                    worst_scores[seq_id] = score
    except Exception as e:
        print(f"Erreur lors de l'extraction des scores du fichier {aln_file}: {e}")

    return worst_scores

def save_worst_scores(worst_scores, output_file):
    """Enregistre les pires scores dans un fichier."""
    create_directory_if_not_exists(os.path.dirname(output_file))
    
    try:
        with open(output_file, 'w') as f:
            for seq_id, score in worst_scores.items():
                f.write(f"{seq_id}\t{score}\n")
        print(f"Scores enregistrés dans {output_file}")
    except Exception as e:
        print(f"Erreur lors de la sauvegarde des scores dans {output_file}: {e}")

def create_blast_db_and_search(input_files_folder, db_folder, aln_folder, result_folder, tmp_folder, max_workers=4):
    input_files = os.listdir(input_files_folder)

    # Utilisation de ThreadPoolExecutor pour paralléliser le traitement des fichiers
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(process_file, filename, input_files_folder, db_folder, aln_folder, result_folder, tmp_folder, max_workers)
            for filename in input_files
        ]

        # Attendre que tous les processus soient terminés
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Erreur lors du traitement d'un fichier : {e}")

# Exemple d'utilisation du script avec sys.argv
if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script.py <input_files_folder> <db_folder> <aln_folder> <result_folder> <tmp_folder>")
        sys.exit(1)

    input_files_folder = sys.argv[1]
    db_folder = sys.argv[2]
    aln_folder = sys.argv[3]
    result_folder = sys.argv[4]
    tmp_folder = sys.argv[5]

    create_blast_db_and_search(input_files_folder, db_folder, aln_folder, result_folder, tmp_folder)
