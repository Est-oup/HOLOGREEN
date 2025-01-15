import os
import subprocess

def sort_bam_by_readname(bam_dir):
    """Trie tous les fichiers BAM dans un répertoire par nom de lecture."""
    # Parcourir tous les fichiers BAM dans le répertoire
    for bam_file in os.listdir(bam_dir):
        if bam_file.endswith('.bam'):
            # Chemin complet vers le fichier BAM
            bam_path = os.path.join(bam_dir, bam_file)
            
            # Fichier BAM trié en sortie
            sorted_bam = os.path.join(bam_dir, f"{bam_file[:-4]}_sorted.bam")
            
            # Commande samtools pour trier par nom de lecture
            print(f"Tri du fichier BAM par nom de lecture : {bam_file}")
            command = f"samtools sort -@ 39 -n -o {sorted_bam} {bam_path}"
            
            # Exécuter la commande
            subprocess.run(command, shell=True, check=True)
            print(f"Fichier BAM trié : {sorted_bam}")

if __name__ == "__main__":
    # Répertoire contenant les fichiers BAM
    bam_dir = "out/bam"
    
    # Trier tous les fichiers BAM par nom de lecture
    sort_bam_by_readname(bam_dir)
