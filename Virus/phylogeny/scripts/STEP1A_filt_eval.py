import os
import sys
import shutil

def filter_and_extract_ids(alignment_folder, output_folder, evalue_thresholds):
    evalue_thresholds = [float(threshold) for threshold in evalue_thresholds.split(',')]
    os.makedirs(output_folder, exist_ok=True)

    for threshold in evalue_thresholds:
        threshold_folder = os.path.join(output_folder, f"{threshold:.0e}")
        if os.path.exists(threshold_folder):
            shutil.rmtree(threshold_folder)
        os.makedirs(threshold_folder, exist_ok=True)

        for filename in os.listdir(alignment_folder):
            filepath = os.path.join(alignment_folder, filename)
            protein_name = filename

            with open(filepath, 'r') as file:
                lines = file.readlines()

            filtered_ids = set()
            for line in lines:
                parts = line.strip().split('\t')
                
                # Vérification du nombre de colonnes pour déterminer la colonne du contig
                if len(parts) > 10:
                    if filename.endswith('_filtered'):  # Pour les fichiers filtrés
                        orf_id = parts[0]  # Colonne 1
                    else:  # Pour les autres fichiers
                        orf_id = parts[1]  # Colonne 2

                    # Filtrage basé sur le seuil d'e-value
                    if float(parts[10]) <= threshold:
                        filtered_ids.add(orf_id)

            output_file = os.path.join(threshold_folder, f"{protein_name}.txt")
            with open(output_file, 'w') as out_file:
                for orf_id in sorted(filtered_ids):
                    out_file.write(f"{orf_id}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_by_threshold.py <alignment_folder> <output_folder> <evalue_thresholds>")
        sys.exit(1)

    alignment_folder = sys.argv[1]
    output_folder = sys.argv[2]
    evalue_thresholds = sys.argv[3]

    filter_and_extract_ids(alignment_folder, output_folder, evalue_thresholds)
