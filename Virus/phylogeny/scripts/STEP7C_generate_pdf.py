import os
import sys
from ete3 import Tree, TreeStyle, NodeStyle, TextFace

def read_reference_sequences(fasta_directory):
    references = set()
    for fasta_file in os.listdir(fasta_directory):
        if fasta_file.endswith(".prt"):
            with open(os.path.join(fasta_directory, fasta_file), 'r') as file:
                for line in file:
                    if line.startswith('>'):
                        sequence_name = line.split()[0][1:]
                        references.add(sequence_name)
    return references

def layout(node, references):
    if node.is_leaf():
        contig_name = node.name
        color = "red" if contig_name in references else "black"
        name_face = TextFace(contig_name, fgcolor=color)
        node.add_face(name_face, column=0, position="branch-right")

def generate_circular_tree(tree_file, output_folder, references):
    try:
        output_file = os.path.join(output_folder, os.path.basename(tree_file).replace(".newick", "_circular.png"))
        tree = Tree(tree_file, format=1)
        tree.set_outgroup(tree.get_midpoint_outgroup())
        
        # TreeStyle configuration for a simple circular tree
        ts = TreeStyle()
        ts.mode = "c"
        ts.show_leaf_name = True
        ts.layout_fn = lambda n: layout(n, references)
        
        tree.render(output_file, w=800, units="px", tree_style=ts)
        print(f"Circular tree PNG generated for {tree_file}")
        
    except Exception as e:
        print(f"Error generating circular tree for {tree_file}: {e}", file=sys.stderr)

def main(tree_folder, output_folder, fasta_directory):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    references = read_reference_sequences(fasta_directory)
    
    tree_files = []
    for root, dirs, files in os.walk(tree_folder):
        for file in files:
            if file.endswith(".newick"):
                tree_files.append(os.path.join(root, file))

    print(f"Found {len(tree_files)} tree files to process.")
    
    for tree_file in tree_files:
        generate_circular_tree(tree_file, output_folder, references)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <tree_folder> <output_folder> <reference_fasta_directory>")
        sys.exit(1)
    
    tree_folder = sys.argv[1]
    output_folder = sys.argv[2]
    fasta_directory = sys.argv[3]

    print(f"Starting circular tree generation with the following parameters:\nTree folder: {tree_folder}\nOutput folder: {output_folder}\nReference FASTA directory: {fasta_directory}")
    main(tree_folder, output_folder, fasta_directory)
