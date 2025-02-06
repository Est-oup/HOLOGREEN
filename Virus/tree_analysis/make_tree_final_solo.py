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

# Charger les groupes de contigs depuis un fichier
def load_contig_groups(filename):
    contig_groups = {}
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            contig, group = line.strip().split('\t')
            contig_groups[contig] = group
    return contig_groups

# Couleurs pour chaque groupe
group_colors = {
    "Ascoviridae": "#E41A1C",
    "Asfarviridae": "#377EB8",
    "Iridoviridae": "#4DAF4A",
    "Marseilleviridae": "#FF7F00",
    "Mimiviridae": "#FF1493",
    "Schizomimiviridae": "#A65628",
    "Mesomimiviridae": "#F781BF",
    "Allomimiviridae": "#999999",
    "Phycodnaviridae": "#A6CEE3",
    "Pitho_like_viruses": "#FFB300",
    "Pandoraviridae": "#33A02C",
    "Unknown": "#000000",
    "Polinton": "#FF4081",
    "Polintonviridae": "#8B4513",
    "Polinton-like": "#FFD700",
    "Lavidaviridae": "#40E0D0",
    "Outgroup": "#7B68EE",
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
def generate_tree(newick_file, output_file, contig_groups, special_contigs):
    if not os.path.exists(newick_file):
        raise FileNotFoundError(f"Le fichier Newick spécifié n'existe pas: {newick_file}")

    try:
        tree = Tree(newick_file, format=1)
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

    # Style de l'arbre
    ts = TreeStyle()
    ts.mode = "c"
    ts.show_branch_support = False
    ts.min_leaf_separation = 2
    ts.arc_start = -180
    ts.arc_span = 360
    ts.margin_right = 5
    ts.margin_left = 5
    ts.margin_top = 10
    ts.margin_bottom = 10
    ts.branch_vertical_margin = -48

    # Fonction pour les étiquettes radiales avec symboles spéciaux
    def radial_labels_layout(node):
        if node.is_leaf():
            group = contig_groups.get(node.name, "Unknown")
            color = group_colors.get(group, "black")
            label_text = node.name if not node.name[0].isdigit() else ""
            label_size = 35
            name_face = TextFace(label_text, fsize=label_size, fgcolor=color, ftype="Arial", bold=True)
            node.add_face(name_face, column=0, position="branch-right")

            # Ajouter les symboles si présents
            symbols = []
            # if node.name in special_contigs['MCP-NCLDV']:
            #     symbols.append(("⚪", "blue"))
            # if node.name in special_contigs['MCP-PLV']:
            #     symbols.append(("⚪", "blue"))
            # if node.name in special_contigs['MCP-virophage']:
            #     symbols.append(("⚪", "blue"))
            # if node.name in special_contigs['MCP-polinton']:
            #     symbols.append(("⚪", "blue"))
            if node.name in special_contigs['dnapol']:
                symbols.append(("⬤", "red"))
            if node.name in special_contigs['ATPase-NCLDV']:
                symbols.append(("✱", "green"))
            if node.name in special_contigs['ATPase-PLV']:
                symbols.append(("✱", "green"))
            if node.name in special_contigs['ATPase-virophage']:
                symbols.append(("✱", "green"))
            if node.name in special_contigs['ATPase-polinton']:
                symbols.append(("✱", "green"))
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

            for i, (shape, sym_color) in enumerate(symbols):
                symbol_face = TextFace(shape, fsize=25, fgcolor=sym_color, bold=True)
                node.add_face(symbol_face, column=1+i, position="branch-right")

    ts.show_leaf_name = False
    ts.layout_fn = radial_labels_layout
    tree.render(output_file, tree_style=ts, w=2000, h=2000)
    print(f"L'image de l'arbre a été sauvegardée à : {output_file}")


def process_tree(tree, group_file, output, occurrence_table_file):
    occurrence_table = load_occurrence_table(occurrence_table_file)
    special_contigs = identify_special_contigs(occurrence_table)
    contig_groups = load_contig_groups(group_file)
    generate_tree(tree, output, contig_groups, special_contigs)

if __name__ == "__main__":
    tree = r'G:\Virus\NCLDV\viral_marker_annotation\phylogeny\out\STEP_7\tree_data_filt\MCP-NCLDV.contree'
    group_file = r'G:\Virus\NCLDV\arbres_rigoux\tree\tree_group\MCP-NCLDV_contig_groups.txt'
    output = r'G:\Virus\NCLDV\arbres_rigoux\tree\tree_svg\MCP-NCLDV.svg'
    occurrence_table_file = r'G:\Virus\NCLDV\arbres_rigoux\tree\presence_absence_table.csv'

    process_tree(tree, group_file, output, occurrence_table_file)
