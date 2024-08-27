import sys
import os
from Bio import SeqIO
from multiprocessing import Pool, Manager

def load_sequences(query_folder):
    """
    Charge toutes les séquences des fichiers FASTA dans un dictionnaire en mémoire.
    """
    sequences = {}
    for file_sequence in os.listdir(query_folder):
        fasta_file = os.path.join(query_folder, file_sequence)
        if os.path.isfile(fasta_file):
            with open(fasta_file, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    sequences[record.id] = record
    return sequences

def extract_sequences_from_ids(args):
    """
    Extrait les séquences correspondant aux IDs fournis et les écrit dans un fichier de sortie.
    """
    ID_list, sequences, output_file = args
    sequences_to_write = [seq for seq_id, seq in sequences.items() if any(partial_id in seq_id for partial_id in ID_list)]

    with open(output_file, "w") as out:
        SeqIO.write(sequences_to_write, out, "fasta")

    return f"Extraction terminée, séquences écrites dans {output_file}"

def extract_sequences_for_orfs(query_folder, ID_folder, out_seqs_folder, max_workers=None):
    os.makedirs(out_seqs_folder, exist_ok=True)

    # Charger toutes les séquences en mémoire
    sequences = load_sequences(query_folder)
    print(f"Total de séquences chargées: {len(sequences)}")

    tasks = []
    for file_name in os.listdir(ID_folder):
        ID_file = os.path.join(ID_folder, file_name)
        if os.path.isfile(ID_file):
            output_file = os.path.join(out_seqs_folder, f"{os.path.splitext(file_name)[0]}_sequences.fasta")

            # Charger les IDs
            with open(ID_file, "r") as f:
                ID_list = {line.strip() for line in f if line.strip()}

            tasks.append((ID_list, sequences, output_file))

    # Utiliser multiprocessing pour paralléliser l'extraction
    with Pool(processes=max_workers) as pool:
        results = pool.map(extract_sequences_from_ids, tasks)
        for result in results:
            print(result)

if __name__ == "__main__":
    # Vérifiez que le nombre d'arguments est correct
    if len(sys.argv) != 4:
        print("Usage: python extractfastafromID.py <query_folder> <ID_folder> <out_seqs_folder>")
        sys.exit(1)

    # Récupérer les chemins des dossiers en entrée
    query_folder = sys.argv[1]
    ID_folder = sys.argv[2]
    out_seqs_folder = sys.argv[3]

    # Déterminer le nombre de processus disponibles
    num_workers = os.cpu_count()

    # Appeler la fonction pour extraire les séquences correspondantes aux IDs
    extract_sequences_for_orfs(query_folder, ID_folder, out_seqs_folder, max_workers=num_workers)
