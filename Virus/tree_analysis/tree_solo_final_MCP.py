import pandas as pd
import os
from ete3 import Tree, NodeStyle, TreeStyle, TextFace

# Charger la table d'occurrence
def load_occurrence_table(filename):
    df = pd.read_csv(filename, encoding='utf-8')
    return df

# Identifier les contigs spéciaux présents dans des colonnes spécifiques
def identify_special_contigs(df):
    special_contigs = {
        "MCP-NCLDV": df[df['MCP-NCLDV'] == True]['taxa'].tolist(),
        "MCP-polinton": df[df['MCP-polinton'] == True]['taxa'].tolist(),
        "MCP-PLV": df[df['MCP-PLV'] == True]['taxa'].tolist(),
        "MCP-virophage": df[df['MCP-virophage'] == True]['taxa'].tolist(),
        "dnapol": df[df['dnapol'] == True]['taxa'].tolist(),
        "ATPase-NCLDV": df[df['ATPase-NCLDV'] == True]['taxa'].tolist(),
        "ATPase-PLV": df[df['ATPase-PLV'] == True]['taxa'].tolist(),
        "ATPase-polinton": df[df['ATPase-polinton'] == True]['taxa'].tolist(),
        "ATPase-virophage": df[df['ATPase-virophage'] == True]['taxa'].tolist(),
        "primase": df[df['primase'] == True]['taxa'].tolist(),
        "rnapol1": df[df['rnapol1'] == True]['taxa'].tolist(),
        "rnapol2": df[df['rnapol2'] == True]['taxa'].tolist(),
        "tf2s": df[df['tf2s'] == True]['taxa'].tolist(),
        "vltf3": df[df['vltf3'] == True]['taxa'].tolist()
    }
    return special_contigs

# Identifier les séquences de référence
def identify_reference_sequences(tree, occurrence_table):
    all_sequences = set(leaf.name for leaf in tree.iter_leaves())
    occurrence_sequences = set(occurrence_table['taxa'].tolist())
    reference_sequences = all_sequences - occurrence_sequences  # Différence d'ensemble
    return reference_sequences

# Charger les groupes de contigs depuis un fichier
def load_contig_groups(filename):
    contig_groups = {}
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            contig, group = line.strip().split('\t')
            contig_groups[contig] = group
    return contig_groups

group_colors = {
    "Pitho_like_viruses": "#90B4EF",
    "Pandoraviridae": "#FFAE42",
    "Polinton": "#9AEF79",
    "Allomimiviridae": "#90E0EF",
    "Missassigned": "#9D4EDD",
    "Polinton-like": "#FFC6FF",
    "Iridoviridae": "#FFD700",
    "Schizomimiviridae": "#582900",  
    "Lavidaviridae":  "#5CDB95",
    "Mimiviridae": "#F56D91",
    "Mesomimiviridae": "#FFDAC1",
    "Phycodnaviridae": "#D48FDE",
    "Asfarviridae": "#3388BB",
    "Unknown": "grey",
    "Ascoviridae": "grey40",
    "Marseilleviridae" : "#FFA07A",
    "Outgroup": "#7B68EE"  # Pas de couleur correspondante en R, laissé inchangé
}




# Enraciner l'arbre selon l'ordre de priorité
def root_tree(tree):
    # Priorité 1 : "Asfarviridae"
    asfarviridae = [node for node in tree.iter_leaves() if "Asfarviridae" in node.name]
    if asfarviridae:
        tree.set_outgroup(asfarviridae[0])
        return

    # Priorité 2 : "Outgroup"
    outgroup_nodes = [node for node in tree.iter_leaves() if "Outgroup" in node.name]
    if outgroup_nodes:
        tree.set_outgroup(outgroup_nodes[0])
        return

    # Si aucun groupe préféré n'est trouvé, utiliser midpoint rooting
    tree.set_outgroup(tree.get_midpoint_outgroup())

# Fonction pour colorer les branches et ajouter des étiquettes
def color_and_label_node(node, color):
    node_style = NodeStyle()
    node_style["fgcolor"] = color
    node_style["size"] = -0
    node_style["vt_line_color"] = color
    node_style["hz_line_color"] = color
    node_style["vt_line_width"] = 6
    node_style["hz_line_width"] = 6
    node.set_style(node_style)

# Fonction principale pour générer l'arbre coloré
def generate_tree(tree, output_file, contig_groups, special_contigs, reference_sequences):
    # Si l'arbre n'est pas un objet Tree, le charger
    if isinstance(tree, str):
        if not os.path.exists(tree):
            raise FileNotFoundError(f"Le fichier Newick spécifié n'existe pas: {tree}")
        try:
            tree = Tree(tree, format=1)
        except Exception as e:
            raise ValueError(f"Erreur lors du chargement du fichier Newick: {e}")
    
    # Enraciner l'arbre avec priorité
    root_tree(tree)
        
    # Fonction pour attribuer les groupes aux nœuds
    def assign_groups(node):
        if node.is_leaf():
            group = contig_groups.get(node.name, "Unknown")
            node.add_feature("group", group)
        else:
            child_groups = set(child.group for child in node.get_children() if hasattr(child, "group"))
            if len(child_groups) == 1:
                node.add_feature("group", child_groups.pop())
            else:
                node.add_feature("group", "mixed")

    # Parcourir l'arbre pour attribuer les groupes
    for node in tree.traverse("postorder"):
        assign_groups(node)

    # Fonction pour colorer les nœuds en fonction du groupe et ajouter des points de bootstrap
    def color_nodes(node):
        group = getattr(node, "group", "Unknown")
        color = group_colors.get(group, "black")
        if group != "mixed":
            color_and_label_node(node, color)

    for node in tree.traverse():
        color_nodes(node)

    # # Style de l'arbre
    # ts = TreeStyle()
    # ts.mode = "c"
    # ts.show_branch_support = False
    # ts.min_leaf_separation = 2
   # ts.arc_start = -180
    # ts.arc_span = 360
    # ts.margin_right = 5
    # ts.margin_left = 5
    # ts.margin_top = 10
    # ts.margin_bottom = 10
    # ts.branch_vertical_margin = -48 

    # Fonction pour les étiquettes radiales avec symboles spéciaux
    def radial_labels_layout(node):
        if node.is_leaf():
            group = contig_groups.get(node.name, "Unknown")
            color = group_colors.get(group, "black")
            
            # Vérifier si le nœud est une séquence de référence
            if node.name in reference_sequences:
                # Ajuster la taille de la croix
                symbol_face = TextFace("✚", fsize=39, fgcolor="black", bold=True)  # Taille réduite pour la croix
                node.add_face(symbol_face, column=0, position="branch-right")
            else:
                # Texte d'étiquette
                label_text = node.name if not node.name[0].isdigit() else ""
                label_size = 49  # Taille standard pour les étiquettes
                name_face = TextFace(label_text, fsize=label_size, fgcolor=color, ftype="Arial", bold=True)
                node.add_face(name_face, column=0, position="branch-right")
                
                # Ajouter les symboles si présents dans les contigs spéciaux
                symbols = []
                if node.name in special_contigs['dnapol']:
                    symbols.append(("⬤", "red"))
                if node.name in special_contigs['ATPase-NCLDV']:
                    symbols.append(("✱", "green"))
                if node.name in special_contigs['ATPase-PLV']:
                    symbols.append(("✱", "blue"))
                if node.name in special_contigs['ATPase-virophage']:
                    symbols.append(("✱", "orange"))
                if node.name in special_contigs['ATPase-polinton']:
                    symbols.append(("✱", "purple"))
                if node.name in special_contigs['primase']:
                    symbols.append(("◆", "cyan"))
                if node.name in special_contigs['rnapol1']:
                    symbols.append(("▲", "green"))
                if node.name in special_contigs['rnapol2']:
                    symbols.append(("▲", "red"))
                if node.name in special_contigs['tf2s']:
                    symbols.append(("★", "purple"))
                if node.name in special_contigs['vltf3']:
                    symbols.append(("✤", "orange"))
                
                # Ajouter les symboles au nœud avec une taille standard
                for i, (shape, sym_color) in enumerate(symbols):
                    symbol_face = TextFace(shape, fsize=20, fgcolor=sym_color, bold=True)
                    node.add_face(symbol_face, column=1 + i, position="branch-right")

    
    # Style de l'arbre avec espacement ajusté
    ts = TreeStyle()
    ts.mode = "c"  # Circulaire
    ts.show_branch_support = False
    ts.min_leaf_separation = 0  # Augmente l'espacement minimal entre les feuilles
    ts.arc_start = -180
    ts.arc_span = 360
    ts.margin_right = 5  # Augmente les marges
    ts.margin_left = 5
    ts.margin_top = 20
    ts.margin_bottom = 20
    ts.branch_vertical_margin = -65  # Évite que les branches soient trop serrées
    
    ts.show_leaf_name = False
    ts.layout_fn = radial_labels_layout
    tree.render(output_file, tree_style=ts, w=2000, h=2000)
    print(f"L'image de l'arbre a été sauvegardée à : {output_file}")



def process_tree(tree_file, group_file, output, occurrence_table_file):
    # Charger l'arbre en tant qu'objet Tree
    try:
        tree = Tree(tree_file, format=1)
    except Exception as e:
        raise ValueError(f"Erreur lors du chargement du fichier Newick: {e}")

    # Charger la table d'occurrence
    occurrence_table = load_occurrence_table(occurrence_table_file)
    
    # Identifier les contigs spéciaux
    special_contigs = identify_special_contigs(occurrence_table)
    
    # Charger les groupes de contigs
    contig_groups = load_contig_groups(group_file)
    
    # Identifier les séquences de référence
    reference_sequences = identify_reference_sequences(tree, occurrence_table)
    
    # Générer l'arbre coloré
    generate_tree(tree, output, contig_groups, special_contigs, reference_sequences)


if __name__ == "__main__":
    tree = r'G:\Virus\NCLDV\viral_marker_annotation\new_phylogeny_iqtree\out\STEP_7\tree_data_filt\MCP-NCLDV.contree'
    group_file = r'G:\Virus\NCLDV\arbres_rigoux_iqtree\tree_group\MCP-NCLDV_contig_groups.txt'
    output = r'G:\Virus\NCLDV\arbres_rigoux_iqtree\tree_svg\MCP-NCLDV.svg'
    occurrence_table_file = r'G:\Virus\NCLDV\arbres_rigoux_iqtree\presence_absence_table.csv'

    process_tree(tree, group_file, output, occurrence_table_file)
