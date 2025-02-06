import os
import pandas as pd
import sys
from Bio import SeqIO

def process_alignment_file(input_file, output_alignment_file, fasta_file, output_fasta_contig_file, output_fasta_bin_file, orf_column, cluster_mapping_file):
    # Lire le fichier d'alignement dans un DataFrame
    df = pd.read_csv(input_file, sep='\t', header=None)
    df.columns = ["orf", "sequence", "score1", "score2", "score3", "score4", "start1", "end1", "start2", "end2", "evalue", "alignment_score"]
    
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
    
    # Combiner les résultats pour les ORFs filtrés au niveau des contigs
    filtered_df_contig = pd.concat([unique_df, best_scores_df])
    
    # Sauvegarder les résultats d'alignement "par contig"
    filtered_df_contig.to_csv(output_alignment_file, sep='\t', index=False, header=False)

    # Extraire les ORFs sélectionnés "par contig"
    selected_orfs_contig = filtered_df_contig['orf'].tolist()

    # Lire les séquences FASTA et extraire celles correspondant aux ORFs sélectionnés
    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    selected_sequences_contig = [seq for seq in fasta_sequences if seq.id.split()[0].split('#')[0] in selected_orfs_contig]

    # Écrire les séquences filtrées "par contig" dans un fichier FASTA
    SeqIO.write(selected_sequences_contig, output_fasta_contig_file, 'fasta')

    # Charger le fichier de mapping cluster/contig
    cluster_mapping = pd.read_csv(cluster_mapping_file, sep='\t', header=None, names=["cluster", "contig"])
    
    # Ajouter l'information des clusters au DataFrame filtré par contig
    filtered_df_contig = pd.merge(filtered_df_contig, cluster_mapping, on="contig", how="inner")

    # Filtrer pour conserver uniquement le meilleur ORF par cluster/bin
    best_orf_per_cluster = filtered_df_contig.loc[filtered_df_contig.groupby('cluster')['alignment_score'].idxmax()]

    # Extraire les ORFs sélectionnés "par cluster/bin"
    selected_orfs_bin = best_orf_per_cluster['orf'].tolist()

    # Réinitialiser le générateur pour lire les séquences FASTA à nouveau
    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    selected_sequences_bin = [seq for seq in fasta_sequences if seq.id.split()[0].split('#')[0] in selected_orfs_bin]

    # Écrire les séquences filtrées "par cluster/bin" dans un fichier FASTA
    SeqIO.write(selected_sequences_bin, output_fasta_bin_file, 'fasta')

def process_all_files(alignment_dir, fasta_dir, cluster_mapping_file, output_alignment_dir, output_fasta_contig_dir, output_fasta_bin_dir):
    for alignment_file in os.listdir(alignment_dir):
        input_file = os.path.join(alignment_dir, alignment_file)
        base_name = os.path.splitext(alignment_file)[0]
        
        # Déterminer la colonne correcte pour les IDs des ORF
        orf_column = 1 if 'filtered' in alignment_file else 0
        
        # Ajouter l'extension .fasta aux fichiers FASTA attendus
        if base_name.startswith("db_"):
            base_name = base_name.replace("db_", "").split('_vs_')[0]
        
        fasta_file = os.path.join(fasta_dir, f"{base_name}.fasta")
        
        if not os.path.exists(fasta_file):
            print(f"FASTA file not found: {fasta_file}, skipping.")
            continue
        
        output_alignment_file = os.path.join(output_alignment_dir, f"{base_name}_filtered.out")
        output_fasta_contig_file = os.path.join(output_fasta_contig_dir, f"{base_name}.fasta")
        output_fasta_bin_file = os.path.join(output_fasta_bin_dir, f"{base_name}.fasta")
        
        print(f"Processing alignment file: {input_file}, corresponding FASTA file: {fasta_file}")
        process_alignment_file(input_file, output_alignment_file, fasta_file, output_fasta_contig_file, output_fasta_bin_file, orf_column, cluster_mapping_file)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python script.py <alignment_dir> <fasta_dir> <cluster_mapping_file> <output_alignment_dir> <output_fasta_contig_dir> <output_fasta_bin_dir>")
        sys.exit(1)
    
    alignment_dir = sys.argv[1]
    fasta_dir = sys.argv[2]
    cluster_mapping_file = sys.argv[3]
    output_alignment_dir = sys.argv[4]
    output_fasta_contig_dir = sys.argv[5]
    output_fasta_bin_dir = sys.argv[6]
    
    # Créer les répertoires de sortie s'ils n'existent pas
    if not os.path.exists(output_alignment_dir):
        os.makedirs(output_alignment_dir)
    
    if not os.path.exists(output_fasta_contig_dir):
        os.makedirs(output_fasta_contig_dir)

    if not os.path.exists(output_fasta_bin_dir):
        os.makedirs(output_fasta_bin_dir)
    
    # Traiter tous les fichiers
    process_all_files(alignment_dir, fasta_dir, cluster_mapping_file, output_alignment_dir, output_fasta_contig_dir, output_fasta_bin_dir)
