import os
import pandas as pd

# Définir les groupes de séquences de référence
ref_groups = {
"Ascoviridae": ['asco.dipuas', 'asco.hevias', 'asco.spofru', 'asco.tricni'],
"Asfarviridae": ['asfa.afswfe', 'asfa.faustv', 'uncl.kaumoe', 'uncl.pacman'],
"Iridoviridae-a": ['irid.lydich', 'irid.lydi01', 'irid.redsea', 'irid.inspki', 'irid.grouir', 'irid.singro','irid.frogvi', 'irid.amtivi'],
"Iridoviridae-b": ['irid.anmiir', ' irid.irvi03', ' irid.wisiri', 'irid.inir25', 'irid.irvi30', 'irid.inir22', 'irid.irvi06', 'irid.arvuir'],
"Marseilleviridae": ['mars.cannvi', 'mars.lausan', 'mars.marsei', 'mars.melbou', 'mars.pormio', 'mars.tunifo'],
"Mimiviridae-related": ['phyc.orgla1', 'phyc.orgla2', 'phyc.phaglo', 'phyc.chervi'],
"Mimiviridae": ['mimi.carovi', 'uncl.catovi', 'uncl.hokovi', 'uncl.klosne', 'uncl.indivi', 'mimi.tupaso', 'mimi.acpomi',  'mimi.acpole', 'mimi.accama', 'mimi.megcou', 'mimi.megchi', 'mimi.moumon', 'mimi.acpomo'],
"Pandoraviridae": ['uncl.molsib', 'pand.pandul', 'pand.pansal'],
"Pitho_like_viruses": ['uncl.pitsib', 'uncl.cedvir', 'uncl.orpheo'],
"Phycodnaviridae": ['phyc.hetaka', 'phyc.osttau', 'phyc.pabuch', 'phyc.actuch', 'phyc.feldvi', 'phyc.ecsivi', 'phyc.chervi', 'phyc.emhu20', 'phyc.emhu86'],
"Polinton.5": ['Polinton.5.NV'],
"Polinton.2": ['Polinton.2.NV', 'Polinton.2.DR', 'Polinton.2.SP'],
"Polinton.1": ['Polinton.1.DY', 'Polinton.1.DR', 'Polinton.1.TC', 'Polinton.1.SP', 'Polinton.1.HM'],
"Polinton.3": ['Polinton.3.TC'],
"Lavidaviridae": ["Spuntnik.YP.002122381.1", "Sputnik.2.YP.009021067.1", "Sputnik.3.YP.009021087.1", "Zamilon.YP.008859635.1", "Yellowstone.Lake.5.YP.009177804.1", "Yellowstone.Lake.6.YP.009177819.1", "Yellowstone.Lake.7.YP.009177692.1"],
"Polinton-like": ["Phaeocystis.globosa.PLV.Gezel.14T.UYE94489.1", "Pleurochrysis.sp.PLV.AUL80826.1", "Pleurochrysis.sp.PLV.AUD57342.1", "Pleurochrysis.sp.PLV.AUD57332.1", "Phaeocystis.globosa.PLV.UYE94198.1", "Phaeocistis.globosa.PGVV14T.0012.PLV.UYE94489.1", "Phaeocistis.globosa.PLV.PGVV.509141021", "Tetraselmis.PLV.Tetvi472343124"],
"Outgroup": ["out.asco.dipuas", "out.asfa.afswfe", "out.irid.amtivi", "out.mars.cannvi", "out.mimi.accama", "out.pand.pandul", "out.phyc.actuch"],
}

def load_contig_groups(filename):
    """Charge les groupes de contigs depuis un fichier."""
    contig_groups = {}
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            contig, group = line.strip().split('\t')
            contig_groups[contig] = group
    return contig_groups

def collect_all_contigs(group_dir):
    """Collecte tous les contigs à partir des fichiers de groupe, en excluant les séquences de référence."""
    all_contigs = set()
    group_files = {}
    
    # Combine all reference sequences into a single set
    ref_sequences = set()
    for group in ref_groups.values():
        ref_sequences.update(group)
    
    for group_file in os.listdir(group_dir):
        if group_file.endswith("_contig_groups.txt"):
            group_name = group_file.replace("_contig_groups.txt", "")
            group_file_path = os.path.join(group_dir, group_file)
            contig_groups = load_contig_groups(group_file_path)
            # Filter out reference sequences
            filtered_contigs = {k: v for k, v in contig_groups.items() if k not in ref_sequences}
            group_files[group_name] = filtered_contigs
            all_contigs.update(filtered_contigs.keys())
    
    return all_contigs, group_files


def generate_contig_table(all_contigs, group_files):
    contig_data = {contig: {} for contig in all_contigs}
    taxo_column = {}
    log_data = []

    for contig in all_contigs:
        taxo_value = "unmatch"
        taxo_set = set()

        for group_name, contig_groups in group_files.items():
            if contig in contig_groups:
                contig_data[contig][group_name] = contig_groups[contig]
                if contig_groups[contig] != "Unknown" and contig_groups[contig] != "unmatch":
                    taxo_set.add(contig_groups[contig])
                    taxo_value = contig_groups[contig]
            else:
                contig_data[contig][group_name] = "unmatch"

        # Check for taxonomic contradictions
        if len(taxo_set) > 1:
            log_data.append({
                "Contig": contig,
                "Taxonomies détectées": ", ".join(taxo_set)
            })
        
        # Set the TAXO value based on priority
        if taxo_value == "unmatch":
            for group_name in contig_data[contig]:
                if contig_data[contig][group_name] == "Unknown":
                    taxo_value = "Unknown"
                    break
        
        taxo_column[contig] = taxo_value

    # Add the TAXO column to the contig data
    for contig in contig_data:
        contig_data[contig]["TAXO"] = taxo_column[contig]

    return contig_data, log_data

def save_to_csv(contig_data, output_file, log_data, log_file):
    df = pd.DataFrame(contig_data).T  # Transpose the data for correct format
    df.to_csv(output_file, encoding='utf-8')

    if log_data:
        log_df = pd.DataFrame(log_data)
        log_df.to_csv(log_file, index=False, encoding='utf-8')
        print(f"Log des contradictions taxonomiques sauvegardé sous : {log_file}")
    else:
        print("Aucune contradiction taxonomique détectée.")

def main(group_dir, output_file, log_file):
    all_contigs, group_files = collect_all_contigs(group_dir)
    contig_data, log_data = generate_contig_table(all_contigs, group_files)
    save_to_csv(contig_data, output_file, log_data, log_file)
    print(f"Tableau généré et sauvegardé sous : {output_file}")

if __name__ == "__main__":
    group_dir = "G:/Virus/NCLDV/arbres_rigoux_iqtree/tree_group/"  # Dossier contenant les fichiers de groupes
    output_file = "G:/Virus/NCLDV/arbres_rigoux_iqtree/contig_summary.csv"  # Fichier de sortie
    log_file = "G:/Virus/NCLDV/arbres_rigoux_iqtree/contig_taxo_conflicts_log.csv"  # Fichier de log pour les contradictions
    main(group_dir, output_file, log_file)
