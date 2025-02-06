import sys
import os
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
from functools import partial

def load_sequences(fasta_file):
    """
    Charge toutes les séquences du fichier FASTA dans un dictionnaire en mémoire.
    """
    sequences = {}
    print(f"Chargement des séquences depuis {fasta_file}...")
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences[record.id] = record
    return sequences

def process_id_file(args):
    """
    Traite un fichier d'ID en extrayant les séquences correspondantes et les écrit dans le fichier de sortie.
    """
    ID_file, sequences, output_file = args

    # Charger les IDs depuis le fichier
    with open(ID_file, "r") as f:
        ID_list = {line.strip() for line in f if line.strip()}

    # Extraction rapide : cherche uniquement les correspondances exactes d'ID
    sequences_to_write = [sequences[seq_id] for seq_id in ID_list if seq_id in sequences]

    # Ecrire toutes les séquences correspondantes dans un seul fichier de sortie
    if sequences_to_write:
        with open(output_file, "w") as out:
            SeqIO.write(sequences_to_write, out, "fasta")
        return f"Extraction terminée, {len(sequences_to_write)} séquences écrites dans {output_file}"
    else:
        return f"Aucune séquence trouvée pour les IDs dans {output_file}"

def extract_sequences_for_orfs(fasta_file, ID_folder, out_seqs_folder, max_workers=None):
    os.makedirs(out_seqs_folder, exist_ok=True)

    # Charger toutes les séquences du fichier FASTA en mémoire
    sequences = load_sequences(fasta_file)
    print(f"Total de séquences chargées: {len(sequences)}")

    # Préparer les tâches pour les traiter en parallèle
    tasks = []
    for file_name in os.listdir(ID_folder):
        ID_file = os.path.join(ID_folder, file_name)
        if os.path.isfile(ID_file):
            output_file = os.path.join(out_seqs_folder, f"{os.path.splitext(file_name)[0]}.fasta")
            tasks.append((ID_file, sequences, output_file))

    # Utiliser ThreadPool pour gérer les E/S plus efficacement
    with Pool(processes=max_workers) as pool:
        results = pool.map(process_id_file, tasks)
        for result in results:
            print(result)

if __name__ == "__main__":
    # Vérifier que le nombre d'arguments est correct
    if len(sys.argv) != 4:
        print("Usage: python extractfastafromID.py <fichier_fasta> <dossier_ID> <dossier_sortie_sequences>")
        sys.exit(1)

    # Récupérer les chemins des fichiers et dossiers en entrée
    fasta_file = sys.argv[1]
    ID_folder = sys.argv[2]
    out_seqs_folder = sys.argv[3]

    # Déterminer le nombre de processus disponibles
    num_workers = cpu_count()

    # Appeler la fonction pour extraire les séquences correspondantes aux IDs
    extract_sequences_for_orfs(fasta_file, ID_folder, out_seqs_folder, max_workers=num_workers)
