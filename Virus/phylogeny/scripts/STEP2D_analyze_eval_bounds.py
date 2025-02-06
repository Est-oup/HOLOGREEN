import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import re

def count_orfs_and_contigs_in_file(file_path):
    orf_count = 0
    contigs = set()
    
    with open(file_path, 'r') as file:
        for line in file:
            orf_count += 1
            contig_name = line.split('_')[0]
            contigs.add(contig_name)
    
    contig_count = len(contigs)
    return orf_count, contig_count

def count_sequences_in_fasta(fasta_path):
    sequence_count = 0
    
    with open(fasta_path, 'r') as file:
        for line in file:
            if line.startswith(">"):
                sequence_count += 1
    
    print(f"Fichier {fasta_path} contient {sequence_count} séquences.")
    return sequence_count

def process_subdirectory(subdir_path, fasta_dir1, fasta_dir2):
    orf_counts = {}
    contig_counts = {}
    fasta_counts1 = {}
    fasta_counts2 = {}
    
    for file_name in os.listdir(subdir_path):
        if file_name.endswith(".txt"):
            gene_marker = file_name.replace(".txt", "")
            file_path = os.path.join(subdir_path, file_name)
            orf_count, contig_count = count_orfs_and_contigs_in_file(file_path)
            orf_counts[gene_marker] = orf_count
            contig_counts[gene_marker] = contig_count
            
            # Traitement des fichiers FASTA
            fasta_file1 = os.path.join(fasta_dir1, gene_marker + ".fasta")
            fasta_file2 = os.path.join(fasta_dir2, gene_marker + ".fasta")
            
            if os.path.exists(fasta_file1):
                fasta_counts1[gene_marker] = count_sequences_in_fasta(fasta_file1)
            else:
                fasta_counts1[gene_marker] = 0
                print(f"Fichier {fasta_file1} non trouvé.")
                
            if os.path.exists(fasta_file2):
                fasta_counts2[gene_marker] = count_sequences_in_fasta(fasta_file2)
            else:
                fasta_counts2[gene_marker] = 0
                print(f"Fichier {fasta_file2} non trouvé.")
    
    return orf_counts, contig_counts, fasta_counts1, fasta_counts2

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text for text in re.split(r'([+-]?\d+)', s)]

def plot_bars_inverted_sorted(dataframe, title, output_file):
    # Séparer les colonnes numériques (comme "1e+00") des autres
    numeric_cols = [col for col in dataframe.columns if col.startswith('1e')]
    other_cols = [col for col in dataframe.columns if not col.startswith('1e')]

    # Trier les colonnes numériques naturellement
    sorted_numeric_cols = sorted(numeric_cols, key=natural_sort_key)
    
    # Combiner les colonnes triées avec les autres colonnes
    sorted_columns = sorted_numeric_cols + other_cols
    dataframe = dataframe[sorted_columns]
    
    ax = dataframe.T.plot(kind='bar', figsize=(12, 8), width=0.8)
    ax.set_title(title)
    ax.set_xlabel('Gène Marqueur')
    ax.set_ylabel('Nombre d\'occurrences')
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def main(input_directory, output_path, fasta_dir1, fasta_dir2):
    orf_summary = {}
    contig_summary = {}
    fasta_summary1_orf = {}
    fasta_summary2_orf = {}
    fasta_summary1_contig = {}
    fasta_summary2_contig = {}
    
    for subdir_name in os.listdir(input_directory):
        subdir_path = os.path.join(input_directory, subdir_name)
        if os.path.isdir(subdir_path):
            orf_counts, contig_counts, fasta_counts1, fasta_counts2 = process_subdirectory(subdir_path, fasta_dir1, fasta_dir2)
            orf_summary[subdir_name] = orf_counts
            contig_summary[subdir_name] = contig_counts
            fasta_summary1_orf[subdir_name] = fasta_counts1
            fasta_summary2_orf[subdir_name] = fasta_counts2
            fasta_summary1_contig[subdir_name] = fasta_counts1
            fasta_summary2_contig[subdir_name] = fasta_counts2
    
    # Combine ORF and contig summaries with FASTA counts into a single DataFrame
    combined_orf_df = pd.concat([pd.DataFrame(orf_summary).fillna(0).astype(int).T,
                                 pd.DataFrame(fasta_summary1_orf).fillna(0).astype(int).T,
                                 pd.DataFrame(fasta_summary2_orf).fillna(0).astype(int).T],
                                axis=1, keys=['ORF Count', 'FASTA Count Dir1', 'FASTA Count Dir2'])

    combined_contig_df = pd.concat([pd.DataFrame(contig_summary).fillna(0).astype(int).T,
                                    pd.DataFrame(fasta_summary1_contig).fillna(0).astype(int).T,
                                    pd.DataFrame(fasta_summary2_contig).fillna(0).astype(int).T],
                                   axis=1, keys=['Contig Count', 'FASTA Count Dir1', 'FASTA Count Dir2'])

    # Sort the DataFrame by the index (i.e., subdir_name)
    combined_orf_df = combined_orf_df.sort_index()
    combined_contig_df = combined_contig_df.sort_index()

    # Save the summaries to CSV files
    if os.path.isdir(output_path):
        orf_output_file = os.path.join(output_path, 'orf_summary.csv')
        contig_output_file = os.path.join(output_path, 'contig_summary.csv')
    else:
        orf_output_file = output_path.replace('.csv', '_orf_summary.csv')
        contig_output_file = output_path.replace('.csv', '_contig_summary.csv')
    
    combined_orf_df.to_csv(orf_output_file, index=True)
    combined_contig_df.to_csv(contig_output_file, index=True)
    
    print(f"Les tableaux récapitulatifs ont été sauvegardés dans {orf_output_file} et {contig_output_file}")
    
    # Generate bar plots
    orf_plot_file = orf_output_file.replace('.csv', '.png')
    contig_plot_file = contig_output_file.replace('.csv', '.png')
    
    plot_bars_inverted_sorted(combined_orf_df['ORF Count'], 'Nombre d\'ORF par Gène Marqueur', orf_plot_file)
    plot_bars_inverted_sorted(combined_contig_df['Contig Count'], 'Nombre de Contigs par Gène Marqueur', contig_plot_file)
    
    print(f"Les diagrammes en barres ont été sauvegardés dans {orf_plot_file} et {contig_plot_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <dossier_filt_results> <fichier_ou_dossier_sortie> <dossier_fasta1> <dossier_fasta2>")
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
