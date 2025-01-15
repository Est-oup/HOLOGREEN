import os
from Bio import SeqIO
import subprocess

def extract_bam_order(bam_dir):
    """Extrait l'ordre des contigs du premier fichier BAM trouvé dans le répertoire."""
    bam_files = [f for f in os.listdir(bam_dir) if f.endswith('.bam')]

    if not bam_files:
        raise FileNotFoundError(f"Aucun fichier BAM trouvé dans le répertoire {bam_dir}")
    
    # Prendre le premier fichier BAM
    bam_file = os.path.join(bam_dir, bam_files[0])

    # Extraire l'ordre des contigs avec samtools view
    print(f"Extraction de l'ordre des contigs à partir de {bam_file}")
    command = f"samtools view -H {bam_file} | grep '@SQ' | cut -f 2 | cut -d ':' -f 2"
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, text=True)

    contig_order = result.stdout.strip().split('\n')

    # Sauvegarder l'ordre des contigs dans un fichier texte
    order_file = os.path.join(bam_dir, "bam_contig_order.txt")
    with open(order_file, 'w') as f:
        for contig in contig_order:
            f.write(f"{contig}\n")
    
    print(f"Ordre des contigs sauvegardé dans {order_file}")
    return order_file

def sort_fasta_by_bam_order(fasta_file, order_file, output_fasta):
    """Trie un fichier FASTA selon l'ordre des contigs dans un fichier d'ordre."""
    with open(order_file, 'r') as order_handle:
        contig_order = [line.strip() for line in order_handle]

    fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # Écrire les séquences triées dans un fichier FASTA
    with open(output_fasta, 'w') as output_handle:
        for contig_name in contig_order:
            if contig_name in fasta_sequences:
                SeqIO.write(fasta_sequences[contig_name], output_handle, "fasta")
            else:
                print(f"Contig {contig_name} non trouvé dans {fasta_file}")

if __name__ == "__main__":
    # Chemins vers les fichiers et répertoires
    bam_dir = "out/bam"  # Dossier contenant les fichiers BAM
    fasta_file = "Data/Contigs/HOLOGREEN_contigs_cleaned.fasta"
    output_fasta = "Data/Contigs/HOLOGREEN_contigs_sorted.fasta"
    
    # Extraire l'ordre des contigs du premier fichier BAM
    order_file = extract_bam_order(bam_dir)

    # Trier le fichier FASTA selon l'ordre des contigs du BAM
    sort_fasta_by_bam_order(fasta_file, order_file, output_fasta)
    print(f"Le fichier FASTA trié est sauvegardé dans {output_fasta}")
