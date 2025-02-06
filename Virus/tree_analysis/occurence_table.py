import os
from Bio import Phylo
from io import StringIO
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet

# Dossiers d'entrée et de sortie
input_folder ='G:/Virus/NCLDV/viral_marker_annotation/new_phylogeny_iqtree/out/STEP_7/tree_data_filt/' 
output_file = 'G:/Virus/NCLDV/arbres_rigoux_iqtree/presence_absence_table.csv'
upset_plot_file = 'G:/Virus/NCLDV/arbres_rigoux_iqtree/upset_plot.png'

# Groupes de séquences de référence
groupes = {
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

# Ensemble de toutes les séquences de référence
all_reference_sequences = set(seq for sequences in groupes.values() for seq in sequences)

# Fonction pour lire les fichiers Newick/contree et extraire les taxons
def extract_taxa_from_tree(file_path):
    try:
        with open(file_path, 'r') as file:
            newick_str = file.read()
        tree = Phylo.read(StringIO(newick_str), "newick")
        taxa = set()
        for clade in tree.find_clades():
            if clade.name and clade.name not in all_reference_sequences:
                taxa.add(clade.name)
        return taxa
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier {file_path}: {e}")
        return set()

# Lire tous les fichiers .contree et extraire les taxons
data = {}
file_list = os.listdir(input_folder)
contree_files = [f for f in file_list if f.endswith('.contree')  or f.endswith('.newick')]

for file_name in contree_files:
    # Retirer l'extension .contree pour obtenir uniquement le nom de la protéine
    protein_name = file_name.split('_')[0].replace('.contree', '') 
    file_path = os.path.join(input_folder, file_name)
    taxa = extract_taxa_from_tree(file_path)
    if taxa:
        if protein_name not in data:
            data[protein_name] = set()
        data[protein_name].update(taxa)


# Créer la table de présence/absence
all_taxa = sorted(set().union(*data.values()))

rows = []
for taxon in all_taxa:
    row = {'taxa': taxon}
    for protein, taxa in data.items():
        row[protein] = taxon in taxa
    rows.append(row)

# Créer un DataFrame à partir des données
df = pd.DataFrame(rows)

# Enregistrer la table dans un fichier CSV
df.to_csv(output_file, index=False)

# Convertir les colonnes en booléens pour UpSetPlot
for protein in data.keys():
    df[protein] = df[protein].astype(bool)

# Préparer les données pour UpSetPlot
df_upset = df.set_index(list(df.columns[1:]))

# Générer le plot UpSet
upset = UpSet(df_upset, subset_size='count', show_counts=True)
upset.plot()
plt.suptitle('UpSet Plot des Taxa non Références par Protéine')
plt.savefig(upset_plot_file)

print(f"Table de présence/absence pour les non-références créée et sauvegardée dans {output_file}")
print(f"Plot UpSet pour les non-références créé et sauvegardé dans {upset_plot_file}")
