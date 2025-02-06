import os
import sys
from Bio import SeqIO

# Vérifier que les bons arguments sont fournis
if len(sys.argv) != 5:
    print("Usage: python script.py <orf_codants_dir> <seq_references_dir> <outgroup_dir> <output_dir>")
    sys.exit(1)

# Récupérer les chemins depuis les arguments de la ligne de commande
orf_codants_dir = sys.argv[1]
seq_references_dir = sys.argv[2]
outgroup_dir = sys.argv[3]
output_dir = sys.argv[4]

# Assurez-vous que le répertoire de sortie existe
os.makedirs(output_dir, exist_ok=True)

# Fonction pour combiner les séquences ORF, les séquences de référence et les outgroups dans un fichier de sortie
def combine_orfs_references_and_outgroup(protein_name):
    orf_file_path = os.path.join(orf_codants_dir, f"{protein_name}.fasta")
    
    # Gestion du nom sans le suffixe '_filtered'
    protein_name_base = protein_name.replace('_filtered', '')
    seq_ref_file_path = os.path.join(seq_references_dir, f"{protein_name_base}.fasta")
    outgroup_file_path = os.path.join(outgroup_dir, f"{protein_name_base}.fasta")
    output_file_path = os.path.join(output_dir, f"{protein_name}.fasta")

    with open(output_file_path, 'w') as out_f:
        # Écrire les ORF s'ils existent
        if os.path.exists(orf_file_path):
            print(f"Ajout des ORF depuis {orf_file_path}")
            for record in SeqIO.parse(orf_file_path, "fasta"):
                SeqIO.write(record, out_f, "fasta")

        # Écrire les séquences de référence s'il y en a
        if os.path.exists(seq_ref_file_path):
            print(f"Ajout des séquences de référence depuis {seq_ref_file_path}")
            for record in SeqIO.parse(seq_ref_file_path, "fasta"):
                SeqIO.write(record, out_f, "fasta")

        # Écrire les outgroups s'ils existent
        if os.path.exists(outgroup_file_path):
            print(f"Ajout des outgroups depuis {outgroup_file_path}")
            for record in SeqIO.parse(outgroup_file_path, "fasta"):
                SeqIO.write(record, out_f, "fasta")

# Parcourir tous les fichiers ORF dans le répertoire
for filename in os.listdir(orf_codants_dir):
    if filename.endswith(".fasta"):
        protein_name = os.path.splitext(filename)[0]  # Récupérer le nom du fichier sans l'extension
        combine_orfs_references_and_outgroup(protein_name)

print("Fusion des ORF, séquences de référence et outgroups terminée.")
