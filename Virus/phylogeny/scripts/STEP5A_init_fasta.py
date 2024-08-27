import os
import sys
from Bio import SeqIO

# Vérifier que les bons arguments sont fournis
if len(sys.argv) != 4:
    print("Usage: python script.py <orf_codants_dir> <seq_references_dir> <output_dir>")
    sys.exit(1)

# Récupérer les chemins depuis les arguments de la ligne de commande
orf_codants_dir = sys.argv[1]
seq_references_dir = sys.argv[2]
output_dir = sys.argv[3]

# Assurez-vous que le répertoire de sortie existe
os.makedirs(output_dir, exist_ok=True)

# Fonction pour combiner les séquences ORF et les séquences de référence dans un fichier de sortie
def combine_orfs_and_references(protein_name):
    orf_file_path = os.path.join(orf_codants_dir, f"{protein_name}.fasta")
    
    # Supprimer le suffixe '_filtered' s'il existe
    protein_name_base = protein_name.replace('_filtered', '')
    seq_ref_file_path = os.path.join(seq_references_dir, f"{protein_name_base}.prt")
    output_file_path = os.path.join(output_dir, f"{protein_name}.fasta")

    with open(output_file_path, 'w') as out_f:
        # Écrire les ORF
        if os.path.exists(orf_file_path):
            for record in SeqIO.parse(orf_file_path, "fasta"):
                SeqIO.write(record, out_f, "fasta")

        # Écrire les séquences de référence
        if os.path.exists(seq_ref_file_path):
            for record in SeqIO.parse(seq_ref_file_path, "fasta"):
                SeqIO.write(record, out_f, "fasta")

# Parcourir tous les fichiers ORF dans le répertoire
for filename in os.listdir(orf_codants_dir):
    if filename.endswith(".fasta"):
        protein_name = os.path.splitext(filename)[0]
        combine_orfs_and_references(protein_name)

print("Combinaison des ORF et des séquences de référence terminée.")
