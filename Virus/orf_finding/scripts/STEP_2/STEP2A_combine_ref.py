import os
import sys
from Bio import SeqIO

def combine_fasta(input_folder, output_folder):
    # Vérifie que les dossiers existent
    if not os.path.exists(input_folder):
        print(f"Erreur : Le dossier d'entrée {input_folder} n'existe pas.")
        sys.exit(1)
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Fichier de sortie
    output_file = os.path.join(output_folder, "combined_sequences.fasta")

    # Liste pour stocker toutes les séquences combinées
    combined_sequences = []

    # Parcourt chaque fichier dans le dossier d'entrée
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".fasta") or file_name.endswith(".fa"):  # Filtre les fichiers au format FASTA
            file_path = os.path.join(input_folder, file_name)
            print(f"Traitement du fichier : {file_name}")
            
            # Extrait le nom du fichier (sans l'extension) pour l'ajouter au titre de la séquence
            base_name = os.path.splitext(file_name)[0]

            # Parcours de chaque séquence dans le fichier FASTA
            with open(file_path, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    # Ajout du suffixe _X pour identifier la provenance de la séquence
                    record.id = f"{record.id}_{base_name}"
                    record.description = f"{record.description}_{base_name}"
                    combined_sequences.append(record)

    # Écriture des séquences combinées dans le fichier de sortie
    with open(output_file, "w") as output_handle:
        SeqIO.write(combined_sequences, output_handle, "fasta")
    
    print(f"Les séquences combinées ont été écrites dans {output_file}")

if __name__ == "__main__":
    # Vérifie les arguments de la ligne de commande
    if len(sys.argv) != 3:
        print("Usage : python combine_fasta.py <dossier_d'entrée> <dossier_de_sortie>")
        sys.exit(1)

    # Arguments de la ligne de commande
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    # Appelle la fonction principale
    combine_fasta(input_folder, output_folder)
