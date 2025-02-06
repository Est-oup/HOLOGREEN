import os
import csv
from Bio import SeqIO

# Chemin du dossier contenant les fichiers FASTA
input_folder = "G:/Virus/NCLDV/arbres_rigoux/seq_references/"

# Chemin du fichier CSV de sortie
output_csv = "G:/Virus/NCLDV/arbres_rigoux/seq_references/summary_supp_table.csv"

def fasta_to_csv(input_folder, output_csv):
    # Ouvrir le fichier CSV en mode écriture
    with open(output_csv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        
        # Parcourir tous les fichiers dans le dossier
        for filename in os.listdir(input_folder):
            if filename.endswith(".fasta") or filename.endswith(".fa"):  # Filtrer les fichiers FASTA
                fasta_path = os.path.join(input_folder, filename)
                
                # Ajouter une ligne avec le nom du fichier
                csvwriter.writerow([f"Fichier : {filename}"])
                csvwriter.writerow(["Entête", "Séquence"])  # Entêtes des colonnes

                # Lire le fichier FASTA et extraire les entêtes et les séquences
                with open(fasta_path, "r") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        header = record.id
                        sequence = str(record.seq)
                        csvwriter.writerow([header, sequence])

                # Ajouter une ligne vide entre chaque fichier pour lisibilité
                csvwriter.writerow([])

    print(f"CSV généré avec succès : {output_csv}")

if __name__ == "__main__":
    # Appeler la fonction pour générer le CSV
    fasta_to_csv(input_folder, output_csv)

