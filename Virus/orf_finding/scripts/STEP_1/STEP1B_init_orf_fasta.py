import os
import sys
from Bio import SeqIO

def extract_orf_ids(prodigal_folder):
    orf_ids = set()
    for filename in os.listdir(prodigal_folder):
        if filename.endswith('.fasta'):
            for record in SeqIO.parse(os.path.join(prodigal_folder, filename), 'fasta'):
                orf_id = record.id.split('_')[0]
                orf_ids.add(orf_id)
    return orf_ids

def create_fasta_with_orf_ids(query_folder, orf_ids, output_file, log_file):
    sequences_to_write = []
    extracted_count = 0

    for filename in os.listdir(query_folder):
        for record in SeqIO.parse(os.path.join(query_folder, filename), 'fasta'):
            seq_id = record.id  # Utilise l'ID complet de la séquence
            if any(orf_id in seq_id for orf_id in orf_ids):  # Vérifie si un ID partiel est contenu dans l'ID de la séquence
                sequences_to_write.append(record)
                extracted_count += 1

    # Assurez-vous que le chemin du fichier de sortie est valide et non un répertoire
    if os.path.isdir(output_file):
        output_file = os.path.join(output_file, "output.fasta")

    with open(output_file, 'w') as output_handle:
        SeqIO.write(sequences_to_write, output_handle, 'fasta')

    with open(log_file, 'w') as log_handle:
        log_handle.write(f"Nombre d'ORF trouvés: {len(orf_ids)}\n")
        log_handle.write(f"Nombre de séquences extraites: {extracted_count}\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python extract_and_create_fasta.py <prodigal_folder> <query_folder> <output_file> <log_file>")
        sys.exit(1)

    prodigal_folder = sys.argv[1]
    query_folder = sys.argv[2]
    output_file = sys.argv[3]
    log_file = sys.argv[4]

    orf_ids = extract_orf_ids(prodigal_folder)
    print(f"IDs des ORF extraits: {orf_ids}")  # Ajout d'un message de débogage
    create_fasta_with_orf_ids(query_folder, orf_ids, output_file, log_file)
    print(f"Extraction et création du fichier {output_file} terminées. Détails dans le fichier de log {log_file}.")
