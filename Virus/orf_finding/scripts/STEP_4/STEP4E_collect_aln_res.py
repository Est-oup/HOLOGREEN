# -- coding: utf-8 --

import os
import sys
import pandas as pd
from collections import defaultdict

# Fonction pour lire les résultats d'alignement
def read_alignment_results(aln_folder):
    aln_files = [os.path.join(aln_folder, f) for f in os.listdir(aln_folder) if f.endswith(".m8")]
    results = defaultdict(list)

    for aln_file in aln_files:
        df = pd.read_csv(aln_file, sep='\t', header=None)
        df.columns = ['orf_id', 'ref_seq', 'identity', 'alignment_length', 'mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']
        results[aln_file] = df

    return results

# Fonction pour filtrer les meilleures correspondances
def filter_best_matches(df):
    best_matches = df.loc[df.groupby('orf_id')['bit_score'].idxmax()].reset_index(drop=True)
    return best_matches

# Fonction pour gérer les ambiguïtés
def resolve_ambiguities(df):
    resolved = []
    ambiguous = []
    
    grouped = df.groupby('orf_id')
    for orf_id, group in grouped:
        if len(group['evalue'].unique()) == 1:
            resolved.append(group)
        else:
            ambiguous.append(group)

    if resolved:
        resolved_df = pd.concat(resolved).reset_index(drop=True)
    else:
        resolved_df = pd.DataFrame()

    if ambiguous:
        ambiguous_df = pd.concat(ambiguous).reset_index(drop=True)
    else:
        ambiguous_df = pd.DataFrame()

    return resolved_df, ambiguous_df

# Fonction pour écrire les IDs dans des fichiers
def write_ids(output_folder, file_name, ids):
    with open(os.path.join(output_folder, file_name), 'w') as f:
        for orf_id in sorted(ids):
            f.write(f"{orf_id}\n")

# Fonction pour extraire les séquences fasta
def extract_fasta_sequences(fasta_folder, ids, output_folder, filename):
    output_path = os.path.join(output_folder, filename)
    with open(output_path, 'w') as output_file:
        for fasta_file in os.listdir(fasta_folder):
            if fasta_file.endswith('.fasta'):
                with open(os.path.join(fasta_folder, fasta_file), 'r') as f:
                    sequences = f.read().split('>')[1:]
                    sequences_dict = {seq.split('\n')[0].split()[0]: seq for seq in sequences}
                    matching_sequences = [sequences_dict[orf_id] for orf_id in ids if orf_id in sequences_dict]
                    for seq in matching_sequences:
                        output_file.write(f'>{seq}')

# Fonction pour traiter les fichiers d'alignement
def process_alignment_files(results, filtered_output_folder, final_output_folder, orf_folder, fasta_output_folder):
    # Créer les répertoires de sortie si nécessaire
    if not os.path.exists(filtered_output_folder):
        os.makedirs(filtered_output_folder)
    if not os.path.exists(final_output_folder):
        os.makedirs(final_output_folder)
    if not os.path.exists(fasta_output_folder):
        os.makedirs(fasta_output_folder)
    
    log_file_path = os.path.join(filtered_output_folder, "ambiguities.log")
    log_file = open(log_file_path, 'w')
    
    for aln_file, df in results.items():
        best_matches = filter_best_matches(df)
        resolved_df, ambiguous_df = resolve_ambiguities(best_matches)

        # Enregistrer les résultats filtrés dans le dossier filtered_output_folder
        base_name = os.path.basename(aln_file).replace('.m8', '')
        protein_name = base_name.split('_vs_')[1]  # Récupérer le nom de la protéine
        resolved_output_path = os.path.join(filtered_output_folder, f"{base_name}_filtered.txt")
        ambiguous_output_path = os.path.join(filtered_output_folder, f"{base_name}_ambiguous.txt")

        if not resolved_df.empty:
            resolved_df.to_csv(resolved_output_path, sep='\t', index=False)
        if not ambiguous_df.empty:
            ambiguous_df.to_csv(ambiguous_output_path, sep='\t', index=False)
            # Écrire les ambiguïtés dans le fichier log
            log_file.write(f"Ambiguities in {base_name}:\n")
            log_file.write(ambiguous_df.to_string(index=False))
            log_file.write("\n\n")

        # Compter les ORFs Polinto, non-Polinto, et ambigus
        polinto_count = resolved_df[resolved_df['ref_seq'].str.contains('_polinto')].shape[0]
        nopolinto_count = resolved_df[~resolved_df['ref_seq'].str.contains('_polinto')].shape[0]
        ambiguous_count = ambiguous_df.shape[0]

        # Écrire les informations de comptage dans le fichier log
        log_file.write(f"{base_name} - Polinto: {polinto_count}, Non-Polinto: {nopolinto_count}, Ambiguous: {ambiguous_count}\n")

        # Classifier les ORFs et écrire les résultats dans le dossier final_output_folder
        polinto_ids = resolved_df[resolved_df['ref_seq'].str.contains('_polinto')]['orf_id'].tolist()
        nopolinto_ids = resolved_df[~resolved_df['ref_seq'].str.contains('_polinto')]['orf_id'].tolist()

        aln_output_folder = os.path.join(final_output_folder, base_name)
        if not os.path.exists(aln_output_folder):
            os.makedirs(aln_output_folder)

        write_ids(aln_output_folder, "polinto_ids.txt", polinto_ids)
        write_ids(aln_output_folder, "nopolinto_ids.txt", nopolinto_ids)

        # Extraire les séquences fasta correspondantes avec le nom de la protéine
        extract_fasta_sequences(orf_folder, polinto_ids, fasta_output_folder, f"{protein_name}_polinto.fasta")
        extract_fasta_sequences(orf_folder, nopolinto_ids, fasta_output_folder, f"{protein_name}_nopolinto.fasta")
    
    log_file.close()

# Fonction principale
def main(aln_folder, orf_folder, filtered_output_folder, final_output_folder, fasta_output_folder):
    results = read_alignment_results(aln_folder)
    process_alignment_files(results, filtered_output_folder, final_output_folder, orf_folder, fasta_output_folder)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script.py <alignment_results_folder> <orf_sequences_folder> <filtered_output_folder> <final_output_folder> <fasta_output_folder>")
        sys.exit(1)

    aln_folder = sys.argv[1]
    orf_folder = sys.argv[2]
    filtered_output_folder = sys.argv[3]
    final_output_folder = sys.argv[4]
    fasta_output_folder = sys.argv[5]

    main(aln_folder, orf_folder, filtered_output_folder, final_output_folder, fasta_output_folder)
