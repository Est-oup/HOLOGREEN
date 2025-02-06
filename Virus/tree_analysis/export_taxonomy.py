import os 
import pandas as pd
from Bio import Phylo
from io import StringIO

def load_contig_groups(filename):
    contig_groups = {}
    with open(filename, 'r', encoding='utf-8') as f:  # Spécifier l'encodage UTF-8
        for line in f:
            contig, group = line.strip().split('\t')
            contig_groups[contig] = group
    return contig_groups

def parse_newick(filename):
    contig_groups = {}
    with open(filename, 'r', encoding='utf-8') as f:
        newick_data = f.read()
        tree = Phylo.read(StringIO(newick_data), 'newick')
        for clade in tree.find_clades():
            if clade.name:
                contig_groups[clade.name] = clade.name  # Utiliser le nom du clade comme taxonomie par défaut
    return contig_groups

def collect_all_contigs(group_dir):
    all_contigs = set()
    group_files = {}
    
    for group_file in os.listdir(group_dir):
        group_name = os.path.splitext(group_file)[0]  # Nom sans extension
        group_file_path = os.path.join(group_dir, group_file)
        
        if group_file.endswith("_contig_groups.txt"):
            contig_groups = load_contig_groups(group_file_path)
        elif group_file.endswith(".nwk"):
            contig_groups = parse_newick(group_file_path)
        else:
            continue
        
        group_files[group_name] = contig_groups
        all_contigs.update(contig_groups.keys())
    
    return all_contigs, group_files

def generate_contig_table_and_log(all_contigs, group_files):
    contig_data = {contig: {} for contig in all_contigs}
    taxo_column = {}
    log_entries = []

    for contig in all_contigs:
        taxo_value = "unmatch"  # Default value if no match is found
        taxonomies_found = {}  # Track which taxonomies are found in which files

        for group_name, contig_groups in group_files.items():
            if contig in contig_groups:
                contig_data[contig][group_name] = contig_groups[contig]

                if contig_groups[contig] != "Unknown" and contig_groups[contig] != "unmatch":
                    if contig_groups[contig] not in taxonomies_found:
                        taxonomies_found[contig_groups[contig]] = []
                    taxonomies_found[contig_groups[contig]].append(group_name)

                if contig_groups[contig] != "Unknown" and contig_groups[contig] != "unmatch":
                    taxo_value = contig_groups[contig]
            else:
                contig_data[contig][group_name] = "unmatch"
        
        if len(taxonomies_found) > 1:
            log_entries.append({
                "Contig": contig,
                "Taxonomies": ", ".join(f"{taxo} (in {', '.join(files)})" for taxo, files in taxonomies_found.items())
            })

        if taxo_value == "unmatch":
            for group_name in contig_data[contig]:
                if contig_data[contig][group_name] == "Unknown":
                    taxo_value = "Unknown"
                    break
        
        taxo_column[contig] = taxo_value

    for contig in contig_data:
        contig_data[contig]["TAXO"] = taxo_column[contig]

    return contig_data, log_entries

def save_to_csv(contig_data, log_entries, output_file, log_file):
    df = pd.DataFrame(contig_data).T
    df.to_csv(output_file, encoding='utf-8')

    log_df = pd.DataFrame(log_entries)
    log_df.to_csv(log_file, index=False, encoding='utf-8')

def main(group_dir, output_file, log_file):
    all_contigs, group_files = collect_all_contigs(group_dir)
    contig_data, log_entries = generate_contig_table_and_log(all_contigs, group_files)
    save_to_csv(contig_data, log_entries, output_file, log_file)
    print(f"Tableau généré et sauvegardé sous : {output_file}")
    print(f"Log des contradictions sauvegardé sous : {log_file}")

if __name__ == "__main__":
    group_dir = "G:/Virus/NCLDV/arbres_rigoux_iqtree/tree_group"
    output_file = "G:/Virus/NCLDV/arbres_rigoux_iqtree/contig_summary.csv"
    log_file = "G:/Virus/NCLDV/arbres_rigoux_iqtree/contig_taxo_log.csv"
    main(group_dir, output_file, log_file)