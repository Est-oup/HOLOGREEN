import sys
import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_file(file_path, output_folder):
    try:
        print(f"Processing file: {file_path}")  # Ajout de debug
        # Charger le fichier d'alignement
        align_df = pd.read_csv(file_path, sep='\t', header=None)
        
        # Filtrer les lignes correspondant à "Virus" dans la deuxième colonne
        virus_alignments = align_df[align_df[1] == "Virus"]
        
        print(f"Found {len(virus_alignments)} virus alignments")  # Ajout de debug

        # Supprimer les doublons des ID de première colonne
        virus_alignments = virus_alignments.drop_duplicates(subset=[0])
        
        # Sauvegarder le tableau filtré dans le dossier de sortie
        filtered_file = os.path.join(output_folder, os.path.basename(file_path))
        virus_alignments.to_csv(filtered_file, sep='\t', index=False, header=None)
        
        # Retourner les ID des ORFs filtrés
        return virus_alignments[0].tolist()
    except Exception as e:
        print(f"Erreur lors du traitement du fichier {file_path}: {e}")
        return []

def filter_alignments(input_folder, output_folder, output_folder_ID, max_workers=8):
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(output_folder_ID, exist_ok=True)

    filtered_orfs = set()
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_file, os.path.join(input_folder, file_name), output_folder) 
                   for file_name in os.listdir(input_folder) if os.path.isfile(os.path.join(input_folder, file_name))]
        
        for future in as_completed(futures):
            result = future.result()
            print(f"Filtered ORFs: {result}")  # Ajout de debug
            filtered_orfs.update(result)
    
    filtered_orfs_file = os.path.join(output_folder_ID, "dna_pol_filtered_orfs.txt")
    with open(filtered_orfs_file, 'w') as f:
        for orf_id in sorted(filtered_orfs):
            f.write(orf_id + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python STEP3E_DNA_pol_filt.py input_folder output_folder output_folder_ID")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    output_folder_ID = sys.argv[3]
    
    num_workers = os.cpu_count()
    
    filter_alignments(input_folder, output_folder, output_folder_ID, max_workers=num_workers)
