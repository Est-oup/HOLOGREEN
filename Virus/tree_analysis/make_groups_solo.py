import os
from collections import defaultdict
from Bio import Phylo
from multiprocessing import cpu_count
from joblib import Parallel, delayed
from tqdm import tqdm

# Charger les groupes de séquences de référence
groups = {
    "Ascoviridae": ['asco.dipuas', 'asco.hevias', 'asco.spofru', 'asco.tricni'],
    "Asfarviridae": ['asfa.afswfe', 'asfa.faustv', 'uncl.kaumoe', 'uncl.pacman'],
    "Iridoviridae": ['irid.lydich', 'irid.lydi01', 'irid.redsea', 'irid.inspki', 'irid.grouir', 'irid.singro', 'irid.frogvi', 'irid.amtivi',
                       'irid.anmiir', 'irid.irvi03', 'irid.wisiri', 'irid.inir25', 'irid.irvi30', 'irid.inir22', 'irid.irvi06', 'irid.arvuir'],
    "Marseilleviridae": ['mars.cannvi', 'mars.lausan', 'mars.marsei', 'mars.melbou', 'mars.pormio', 'mars.tunifo'],
    "Mesomimiviridae": ['phyc.orgla1', 'phyc.orgla2', 'phyc.phaglo', 'phyc.chervi'],
    "Allomimiviridae": ['tetreaselmis.virus.1.GCF.003057795.1.ASM305779v1', 'pyramimonas.GCF.029885655.1.ASM2988565v1'],
    "Schizomimiviridae": ['prymnesium.kappa.virus.RF01', 'aureococcus.anophagefferens.NC.024697'],
    "Mimiviridae": ['mimi.carovi', 'uncl.catovi', 'uncl.hokovi', 'uncl.klosne', 'uncl.indivi', 'mimi.tupaso', 'mimi.acpomi', 
                    'mimi.acpole', 'mimi.accama', 'mimi.megcou', 'mimi.megchi', 'mimi.moumon', 'mimi.acpomo', 
                    'chlorella.virus.XWO1.GCA.031310605.1.ASM3131060v1', 'chlorella.virus.XW01.ULY68607.1.putative'],
    "Pandoraviridae": ['uncl.molsib', 'pand.pandul', 'pand.pansal'],
    "Pitho_like_viruses": ['uncl.pitsib', 'uncl.cedvir', 'uncl.orpheo'],
    "Phycodnaviridae": ['phyc.hetaka', 'phyc.osttau', 'phyc.pabuch', 'phyc.actuch', 'phyc.feldvi', 'phyc.ecsivi', 'phyc.chervi', 'phyc.emhu20', 'phyc.emhu86'],
    "Polinton": ['Polinton.5.NV','Polinton.2.NV', 'Polinton.2.DR', 'Polinton.2.SP','Polinton.1.DY', 'Polinton.1.DR', 
                 'Polinton.1.TC', 'Polinton.1.SP', 'Polinton.1.HM','Polinton.3.TC'],
    "Lavidaviridae": ["Spuntnik.YP.002122381.1", "Sputnik.2.YP.009021067.1", "Sputnik.3.YP.009021087.1", "Zamilon.YP.008859635.1", 
                      "Yellowstone.Lake.5.YP.009177804.1", "Yellowstone.Lake.6.YP.009177819.1", "Yellowstone.Lake.7.YP.009177692.1",
                      "Spuntnik.YP.002122364.1", "Sputnik.2.YP.009021051.1", "Sputnik.3.YP.009021071.1", "Zamilon.YP.008859647.1", 
                      "Yellowstone.Lake.5.YP.009177784.1", "Yellowstone.Lake.6.YP.009177816.1", "Yellowstone.Lake.YP.009177673.1",
                      "Chlorella.virophage.ULY68431.1", "Chlorella.virophage.ULY68427.1"],
    "Polinton-like": ["Phaeocystis.globosa.PLV.Gezel.14T.UYE94489.1", "Pleurochrysis.sp.PLV.AUL80826.1", "Pleurochrysis.sp.PLV.AUD57342.1", 
                      "Pleurochrysis.sp.PLV.AUD57332.1", "Phaeocystis.globosa.PLV.UYE94198.1", "Phaeocistis.globosa.PGVV14T.0012.PLV.UYE94489.1", 
                      "Phaeocistis.globosa.PLV.PGVV.509141021", "Tetraselmis.PLV.Tetvi472343124"],
    "Outgroup": ["out.asco.dipuas", "out.asfa.afswfe", "out.irid.amtivi", "out.mars.cannvi", "out.mimi.accama", "out.pand.pandul", "out.phyc.actuch"],
}



# Charger le fichier Newick
def load_newick_file(filepath):
    return Phylo.read(filepath, 'newick')

# Obtenir les séquences de contigs
def get_contigs(tree):
    return [clade.name for clade in tree.find_clades() if clade.name]

# Trouver les clades pour chaque groupe de séquences de référence
def get_clades_for_groups(tree, groups):
    group_clades = {}
    for group, refs in groups.items():
        for clade in tree.find_clades():
            clade_names = {c.name for c in clade.get_terminals()}
            if all(ref in clade_names for ref in refs):
                group_clades[group] = clade
                break
    return group_clades

# Pré-calculer les distances entre tous les contigs et les références
def precalculate_distances(tree, contigs, groups):
    distances = {}
    for contig in contigs:
        distances[contig] = {}
        for group, refs in groups.items():
            for ref in refs:
                if ref not in distances[contig]:
                    try:
                        distance = tree.distance(contig, ref)
                        distances[contig][ref] = distance
                    except ValueError:
                        distances[contig][ref] = float('inf')
    return distances

# Fonction pour assigner un contig à un groupe en utilisant les distances pré-calculées
def assign_contig_to_group(contig, groups, distances):
    best_group = "Unknown"
    min_distance = float('inf')
    for group, refs in groups.items():
        group_distance = min(distances[contig].get(ref, float('inf')) for ref in refs)
        if group_distance < min_distance:
            min_distance = group_distance
            best_group = group
    return contig, best_group if min_distance < 1.0 else "Unknown"

# Assigner les contigs aux groupes en utilisant les distances pré-calculées
def assign_groups_by_ancestor(tree, groups, group_clades):
    contigs = get_contigs(tree)
    distances = precalculate_distances(tree, contigs, groups)
    results = Parallel(n_jobs=-1, backend="loky")(
        delayed(assign_contig_to_group)(contig, groups, distances) for contig in tqdm(contigs, desc="Assigning contigs")
    )
    return dict(results)

# Fonction principale
def main(newick_file, output_file):
    tree = load_newick_file(newick_file)
    group_clades = get_clades_for_groups(tree, groups)
    contig_groups = assign_groups_by_ancestor(tree, groups, group_clades)
    
    # Écriture des résultats dans un fichier
    with open(output_file, 'w') as f:
        for contig, group in contig_groups.items():
            f.write(f"{contig}\t{group}\n")

# Utilisation
if __name__ == "__main__":
    newick_file =  r'G:\Virus\NCLDV\viral_marker_annotation\\phylogeny\\out\\STEP_7\\tree_data_filt\\tf2s.contree'
    output_file = 'G:/Virus/NCLDV/arbres_rigoux/tree/tree_group/tf2s_contig_groups.txt'
    main(newick_file, output_file)