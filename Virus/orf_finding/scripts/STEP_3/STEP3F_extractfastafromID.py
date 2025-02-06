import sys
import os
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed

def extract_sequences_from_file(fasta_file, ID_file, output_file):
    # Nom du fichier ID (sans extension)
    base_name = os.path.splitext(os.path.basename(ID_file))[0]

    # Création d'une liste d'ID à partir du fichier donné
    ID_list = set()
    with open(ID_file) as ID:
        for line in ID:
            line = line.strip()
            if line != "":
                ID_list.add(line)

    # Extraction des séquences correspondantes
    sequences_to_write = []
    with open(fasta_file, "r") as f:
        fasta_sequences = SeqIO.parse(f, "fasta")
        
        # Boucle d'extraction des séquences qui correspondent aux IDs
        for seq in fasta_sequences:
            if seq.id in ID_list:
                sequences_to_write.append(seq)

    # Écriture des séquences extraites dans le fichier de sortie
    with open(output_file, "w") as out:
        SeqIO.write(sequences_to_write, out, "fasta")

    print(f"Extraction terminée pour {ID_file}, séquences écrites dans {output_file}")

def extract_sequences_for_orfs(fasta_file, ID_folder, out_seqs_folder, max_workers=None):
    os.makedirs(out_seqs_folder, exist_ok=True)
    
    tasks = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for file_name in os.listdir(ID_folder):
            ID_file = os.path.join(ID_folder, file_name)
            
            # Vérifier si ID_file est bien un fichier
            if os.path.isfile(ID_file):
                output_file = os.path.join(out_seqs_folder, f"{os.path.splitext(file_name)[0]}.fasta")
                tasks.append(executor.submit(extract_sequences_from_file, fasta_file, ID_file, output_file))
        
        for future in as_completed(tasks):
            future.result()  # Pour gérer les exceptions s'il y en a

if __name__ == "__main__":
    # Vérifiez que le nombre d'arguments est correct
    if len(sys.argv) != 4:
        print("Usage: python extractfastafromID.py <fasta_file> <ID_folder> <out_seqs_folder>")
        sys.exit(1)

    # Récupérer les chemins des fichiers en entrée
    fasta_file = sys.argv[1]
    ID_folder = sys.argv[2]
    out_seqs_folder = sys.argv[3]

    # Déterminer le nombre de threads disponibles
    num_workers = os.cpu_count()

    # Appeler la fonction pour extraire les séquences correspondantes aux IDs
    extract_sequences_for_orfs(fasta_file, ID_folder, out_seqs_folder, max_workers=num_workers)
