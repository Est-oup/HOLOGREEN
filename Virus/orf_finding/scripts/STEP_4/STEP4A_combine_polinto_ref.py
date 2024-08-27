# -*- coding: utf-8 -*-
import sys
import os

def combine_sequences(ref_dir, output_dir):
    # Obtenir la liste des fichiers dans le répertoire de référence
    files = os.listdir(ref_dir)
    
    # Filtrer les fichiers qui ont une déclinaison "polinton"
    polinton_files = [f for f in files if '-polinto' in f]
    
    for polinton_file in polinton_files:
        # Trouver le fichier homologue non-polinton
        non_polinton_file = polinton_file.replace('-polinto.prt', '.prt')
        
        # Créer le chemin complet des fichiers
        polinton_path = os.path.join(ref_dir, polinton_file)
        non_polinton_path = os.path.join(ref_dir, non_polinton_file)
        
        # Définir le fichier de sortie
        protein_name = polinton_file.split('-polinto')[0]
        output_file = os.path.join(output_dir, f"{protein_name}")
        
        with open(output_file, 'w', encoding='utf-8') as outfile:
            # Ajouter les séquences polinton
            if os.path.exists(polinton_path):
                with open(polinton_path, 'r', encoding='utf-8') as infile:
                    for line in infile:
                        if line.startswith(">"):
                            outfile.write(f">{line[1:].strip()}_polinto\n")
                        else:
                            outfile.write(line)
            
            # Ajouter les séquences non-polinton si le fichier existe
            if os.path.exists(non_polinton_path):
                with open(non_polinton_path, 'r', encoding='utf-8') as infile:
                    for line in infile:
                        if line.startswith(">"):
                            outfile.write(f">{line[1:].strip()}_nopolinto\n")
                        else:
                            outfile.write(line)
            else:
                print(f"Warning: Nopolinto file {non_polinton_file} does not exist and will be skipped.")

if __name__ == "__main__":
    ref_dir = sys.argv[1]
    output_dir = sys.argv[2]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    combine_sequences(ref_dir, output_dir)
