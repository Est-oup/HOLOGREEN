import sys
import os

def extract_sequences_for_orfs(alignment_folder, info_folder):
    for file_name in os.listdir(alignment_folder):
        alignment_file = os.path.join(alignment_folder, file_name)
        if os.path.isfile(alignment_file):  
            base_name = os.path.splitext(file_name)[0]

            # Créer un dossier de sortie s'il n'existe pas
            output_file = os.path.join(info_folder, "raw", f"{base_name}_ids.txt")
            output_file_unique = os.path.join(info_folder, "dedup", f"{base_name}_id.txt")

            # Liste pour stocker les IDs
            ids = set()

            # Extraire les ID des fichiers d'alignement
            with open(alignment_file, 'r') as align_file:
                for line in align_file:
                    columns = line.strip().split('\t')
                    if len(columns) > 1:  
                        seq_id = columns[1]
                        ids.add(seq_id)

            # Écrire tous les IDs dans le fichier de sortie
            with open(output_file, 'w') as out_file:
                for seq_id in ids:
                    out_file.write(seq_id + '\n')

            # Écrire les IDs uniques dans le fichier de sortie unique
            with open(output_file_unique, 'w') as out_file_unique:
                for seq_id_unique in sorted(ids):
                    out_file_unique.write(seq_id_unique + '\n')


if __name__ == "__main__":
    # Récupérer les chemins des dossiers en entrée
    alignment_folder = sys.argv[1]
    info_folder = sys.argv[2]

    # Appeler la fonction pour extraire les IDs des fichiers d'alignement
    extract_sequences_for_orfs(alignment_folder, info_folder)

