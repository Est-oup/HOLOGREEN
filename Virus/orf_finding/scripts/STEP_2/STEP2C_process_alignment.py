import os
import sys
import pandas as pd
from collections import defaultdict

def filter_best_alignments(alignment_file, output_filtered_file):
    # Lire le fichier d'alignement
    df = pd.read_csv(alignment_file, sep='\t', header=None)
    
    # Renommer les colonnes pour une meilleure lisibilité
    df.columns = ['ref_seq', 'orf', 'identity', 'alignment_length', 'mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']
    
    # Filtrer pour chaque ORF (colonne 'orf') sur la ligne avec le meilleur score d'alignement (bit_score)
    best_alignments = df.loc[df.groupby('orf')['bit_score'].idxmax()].reset_index(drop=True)
    
    # Écrire les résultats dans le fichier filtré
    best_alignments.to_csv(output_filtered_file, sep='\t', index=False)
    return best_alignments

def split_by_markername(filtered_df, output_dir):
    # Créer le dossier de sortie si nécessaire
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Extraire le markername de la colonne ref_seq
    filtered_df['markername'] = filtered_df['ref_seq'].apply(lambda x: x.split('_')[-1])

    # Filtrer et sauvegarder chaque fichier par markername
    for markername, group in filtered_df.groupby('markername'):
        output_file = os.path.join(output_dir, f"{markername}_alignments.txt")
        group.to_csv(output_file, sep='\t', index=False)
        print(f"Fichier créé : {output_file}")

def extract_contig_info(input_dir, output_dir):
    # Créer le dossier de sortie si nécessaire
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file_name in os.listdir(input_dir):
        file_path = os.path.join(input_dir, file_name)
        df = pd.read_csv(file_path, sep='\t', header=0)
        
        # Extraire les contigs à partir du nom d'ORF
        df['contig'] = df['orf'].apply(lambda x: x.split('_')[0])
        
        # Compter les duplications des contigs
        contig_counts = df['contig'].value_counts().reset_index()
        contig_counts.columns = ['contig', 'duplication_count']

        # Enregistrer le résultat
        output_file = os.path.join(output_dir, f"{file_name.replace('_alignments.txt', '')}_contigs.txt")
        contig_counts.to_csv(output_file, sep='\t', index=False)
        print(f"Fichier de contigs créé : {output_file}")

def extract_orf_names(input_dir, output_dir):
    # Créer le dossier de sortie si nécessaire
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file_name in os.listdir(input_dir):
        file_path = os.path.join(input_dir, file_name)
        df = pd.read_csv(file_path, sep='\t', header=0)
        
        # Extraire les noms d'ORFs (colonne 'orf')
        orf_names = df['orf'].unique()

        # Enregistrer le fichier avec les noms d'ORFs
        output_file = os.path.join(output_dir, f"{file_name.replace('_alignments.txt', '')}")
        with open(output_file, 'w') as f:
            for orf in orf_names:
                f.write(orf + '\n')
        print(f"Fichier de noms d'ORFs créé : {output_file}")

def main():
    if len(sys.argv) != 6:
        print("Usage: python process_alignment.py <input_alignment_folder> <filtered_output_folder> <markername_output_folder> <contig_output_folder> <orf_output_folder>")
        sys.exit(1)

    # Arguments
    input_alignment_folder = sys.argv[1]
    filtered_output_folder = sys.argv[2]
    markername_output_folder = sys.argv[3]
    contig_output_folder = sys.argv[4]
    orf_output_folder = sys.argv[5]

    # Fichier d'alignement à traiter (on prend le premier fichier du dossier en entrée)
    alignment_file = os.path.join(input_alignment_folder, os.listdir(input_alignment_folder)[0])

    # Créer le dossier de sortie pour le fichier filtré si nécessaire
    if not os.path.exists(filtered_output_folder):
        os.makedirs(filtered_output_folder)
    
    # Fichier de sortie filtré
    output_filtered_file = os.path.join(filtered_output_folder, 'filtered_alignment.txt')

    # Étape 1: Filtrer l'alignement en fonction du meilleur bit score
    filtered_df = filter_best_alignments(alignment_file, output_filtered_file)

    # Étape 2: Créer des fichiers sous-jacents basés sur le markername
    split_by_markername(filtered_df, markername_output_folder)

    # Étape 3: Extraire les contigs et compter les duplications
    extract_contig_info(markername_output_folder, contig_output_folder)

    # Étape 4: Extraire les noms d'ORFs
    extract_orf_names(markername_output_folder, orf_output_folder)

if __name__ == "__main__":
    main()
