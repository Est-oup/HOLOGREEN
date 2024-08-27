import os
import pandas as pd
import sys
from Bio import SeqIO

def process_alignment_file(input_file, output_alignment_file, fasta_file, output_fasta_file, orf_column):
    # Lire le fichier d'alignement dans un DataFrame
    df = pd.read_csv(input_file, sep='\t', header=None)
    
    # Déterminer les colonnes en fonction du type de fichier
    if 'filtered' in input_file:
        df.columns = ["orf", "sequence", "score1", "score2", "score3", "score4", "start1", "end1", "start2", "end2", "evalue", "alignment_score"]
    else:
        df.columns = ["sequence", "orf", "score1", "score2", "score3", "score4", "start1", "end1", "start2", "end2", "evalue", "alignment_score"]
    
    # Extraire le nom du contig de la colonne 'orf'
    df['contig'] = df['orf'].apply(lambda x: x.split('_')[0])

    # Trouver les contigs dupliqués et uniques
    contig_counts = df['contig'].value_counts()
    unique_contigs = contig_counts[contig_counts == 1].index
    duplicated_contigs = contig_counts[contig_counts > 1].index
    
    # Conserver les lignes avec des contigs uniques
    unique_df = df[df['contig'].isin(unique_contigs)]
    
    # Conserver les lignes avec les scores les plus élevés pour les contigs dupliqués
    duplicated_df = df[df['contig'].isin(duplicated_contigs)]
    best_scores_df = duplicated_df.loc[duplicated_df.groupby('contig')['alignment_score'].idxmax()]
    
    # Combiner les résultats
    final_df = pd.concat([unique_df, best_scores_df])
    
    # Sauvegarder le résultat dans un nouveau fichier
    final_df.to_csv(output_alignment_file, sep='\t', index=False, header=False)
    
    # Extraire les ORF ids sélectionnés
    selected_orfs = final_df['orf'].tolist()
    
    # Lire les séquences FASTA et extraire celles correspondant aux ORF sélectionnés
    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    selected_sequences = []
    
    for seq in fasta_sequences:
        seq_id = seq.id.split()[0].split('#')[0]  # Extraire l'ID correct sans annotations supplémentaires
        if seq_id in selected_orfs:
            selected_sequences.append(seq)
    
    # Écrire les séquences sélectionnées dans un nouveau fichier FASTA
    SeqIO.write(selected_sequences, output_fasta_file, 'fasta')

def process_all_files(alignment_dir, fasta_dir, output_alignment_dir, output_fasta_dir):
    for alignment_file in os.listdir(alignment_dir):
        input_file = os.path.join(alignment_dir, alignment_file)
        base_name = os.path.splitext(alignment_file)[0]
        
        # Déterminer la colonne correcte pour les IDs des ORF
        orf_column = 1 if 'filtered' in alignment_file else 0
        
        fasta_file = os.path.join(fasta_dir, f"{base_name}.fasta")
        output_alignment_file = os.path.join(output_alignment_dir, f"{base_name}")
        output_fasta_file = os.path.join(output_fasta_dir, f"{base_name}.fasta")
        
        process_alignment_file(input_file, output_alignment_file, fasta_file, output_fasta_file, orf_column)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <alignment_dir> <fasta_dir> <output_alignment_dir> <output_fasta_dir>")
        sys.exit(1)
    
    alignment_dir = sys.argv[1]
    fasta_dir = sys.argv[2]
    output_alignment_dir = sys.argv[3]
    output_fasta_dir = sys.argv[4]
    
    if not os.path.exists(output_alignment_dir):
        os.makedirs(output_alignment_dir)
    
    if not os.path.exists(output_fasta_dir):
        os.makedirs(output_fasta_dir)
    
    process_all_files(alignment_dir, fasta_dir, output_alignment_dir, output_fasta_dir)
