def convert_phylip_to_txt(phy_file, txt_output):
    with open(phy_file, "r") as infile, open(txt_output, "w") as outfile:
        for i, line in enumerate(infile):
            # Écrire chaque ligne du fichier PHYLIP dans le fichier texte
            outfile.write(line)
    
    print(f"Le fichier PHYLIP a été converti en fichier texte et sauvegardé dans {txt_output}")

# Fonction pour extraire les séquences du fichier texte, supprimer les gaps et les enregistrer en format FASTA
def extract_fasta_from_txt(txt_input, fasta_output):
    with open(txt_input, "r") as infile, open(fasta_output, "w") as outfile:
        # Sauter la première ligne (en-tête PHYLIP)
        infile.readline()
        
        for line in infile:
            # Vérifier que la ligne contient une séquence (les lignes doivent avoir deux parties : identifiant et séquence)
            parts = line.split()
            if len(parts) == 2:
                seq_id = parts[0]
                sequence = parts[1]
                # Supprimer les gaps (caractères "-") de la séquence
                cleaned_sequence = sequence.replace("-", "")
                # Écrire au format FASTA
                outfile.write(f">{seq_id}\n{cleaned_sequence}\n")
    
    print(f"Les séquences sans gaps ont été extraites et sauvegardées en format FASTA dans {fasta_output}")

# Fichiers d'entrée et de sortie
phy_file = "F:/Virus/NCLDV/arbres_rigoux/plv_seqs/110_PLV_and_Virophage_MCP-ATPase_alignment.phy"  # Remplacez par le chemin vers votre fichier PHYLIP
txt_output = "F:/Virus/NCLDV/arbres_rigoux/plv_seqs/converted_alignment.txt"  # Fichier de sortie en format texte brut
fasta_output = "F:/Virus/NCLDV/arbres_rigoux/plv_seqs/extracted_sequences_no_gaps.fasta"  # Fichier de sortie au format FASTA

# Convertir PHYLIP en fichier texte
convert_phylip_to_txt(phy_file, txt_output)

# Extraire les séquences du fichier texte, supprimer les gaps et les enregistrer en format FASTA
extract_fasta_from_txt(txt_output, fasta_output)
