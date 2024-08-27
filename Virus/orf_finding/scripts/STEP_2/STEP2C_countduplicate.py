import sys
import os

def count_duplicates(input_folder, output_folder):
    print("Input Folder:", input_folder)
    print("Output Folder:", output_folder)

    for file_name in os.listdir(input_folder):
        input_file = os.path.join(input_folder, file_name)
        if os.path.isfile(input_file):
            print("Processing file:", file_name)
            base_name = os.path.splitext(file_name)[0]

            # Enlever le suffixe "_X" du nom de contig
            contig_counts = {}
            with open(input_file, 'r') as input_f:
                for line in input_f:
                    contig_name = line.strip().split('_')[0]  # Enlever le suffixe "_X"
                    if contig_name in contig_counts:
                        contig_counts[contig_name] += 1
                    else:
                        contig_counts[contig_name] = 1

            output_file = os.path.join(output_folder, f"dupcount_{base_name}.txt")
            print("Output File:", output_file)

            print("Duplicates Counted:", len(contig_counts))
            print("Writing to Output File...")

            # Écrire les résultats dans un fichier
import sys
import os

def count_duplicates(input_folder, output_folder):
    print("Input Folder:", input_folder)
    print("Output Folder:", output_folder)

    for file_name in os.listdir(input_folder):
        input_file = os.path.join(input_folder, file_name)
        if os.path.isfile(input_file):
            print("Processing file:", file_name)
            base_name = os.path.splitext(file_name)[0]

            # Enlever le suffixe "_X" du nom de contig
            contig_counts = {}
            with open(input_file, 'r') as input_f:
                for line in input_f:
                    contig_name = line.strip().split('_')[0]  # Enlever le suffixe "_X"
                    if contig_name in contig_counts:
                        contig_counts[contig_name] += 1
                    else:
                        contig_counts[contig_name] = 1

            output_file = os.path.join(output_folder, f"dupcount_{base_name}.txt")
            print("Output File:", output_file)

            print("Duplicates Counted:", len(contig_counts))
            print("Writing to Output File...")

            # Écrire les résultats dans un fichier
            with open(output_file, 'w') as output_f:
                for contig, count in contig_counts.items():
                    output_f.write(f"{contig}\t{count}\n")

            print("Done with file:", file_name)

if __name__ == "__main__":
    # Récupérer les chemins des dossiers en entrée
    info_folder = sys.argv[1]
    output_folder = sys.argv[2]
    # Appeler la fonction pour compter les contigs dupliqués
    count_duplicates(info_folder, output_folder)

