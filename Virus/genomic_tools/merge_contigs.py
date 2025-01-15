from Bio import SeqIO
from Bio.Seq import Seq

# Charger les contigs et les chevauchements
contigs_file = "bin_19974.fasta"
overlaps_file = "significant_overlaps.csv"

# Charger les contigs
contigs = {record.id: record.seq for record in SeqIO.parse(contigs_file, "fasta")}

# Charger les chevauchements
import pandas as pd
overlaps = pd.read_csv(overlaps_file)

merged_contigs = {}
used_contigs = set()

for _, overlap in overlaps.iterrows():
    query = overlap["query"]
    target = overlap["target"]
    if query in used_contigs or target in used_contigs:
        continue  # Ignorer les contigs déjà fusionnés

    # Fusionner les séquences
    query_seq = contigs[query]
    target_seq = contigs[target]

    # Calculer les positions de fusion
    overlap_length = min(len(query_seq), len(target_seq)) // 2
    merged_seq = query_seq + target_seq[overlap_length:]
    merged_contigs[f"{query}_{target}_merged"] = merged_seq

    used_contigs.update([query, target])

# Ajouter les contigs non utilisés
for contig_id, seq in contigs.items():
    if contig_id not in used_contigs:
        merged_contigs[contig_id] = seq

# Sauvegarder les contigs fusionnés
with open("merged_contigs.fasta", "w") as out_fasta:
    for contig_id, seq in merged_contigs.items():
        out_fasta.write(f">{contig_id}\n{str(seq)}\n")

print("Contigs fusionnés sauvegardés dans 'merged_contigs.fasta'.")
