import numpy as np
import pandas as pd
import sys

# Vérifier que les bons arguments sont fournis
if len(sys.argv) != 3:
    print("Usage: python convert_npz_to_csv.py <input_file.npz> <output_file.csv>")
    sys.exit(1)

# Récupérer les fichiers d'entrée et de sortie depuis les arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Charger le fichier .npz
data = np.load(input_file)

# Convertir les données en DataFrame
df = pd.DataFrame({
    'contig_id': list(data.keys()),
    'bin_contig_size': [data[key] for key in data.keys()]
})

# Sauvegarder en CSV
df.to_csv(output_file, index=False)
