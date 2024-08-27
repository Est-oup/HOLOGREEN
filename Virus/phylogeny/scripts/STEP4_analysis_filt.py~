import os
import sys
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

def count_sequences(fasta_folder):
    count_data = []
    for root, dirs, files in os.walk(fasta_folder):
        for file in files:
            if file.endswith(".fasta"):
                file_path = os.path.join(root, file)
                count = sum(1 for _ in SeqIO.parse(file_path, "fasta"))
                count_data.append({"Protein": os.path.splitext(file)[0], "Count": count})
    return pd.DataFrame(count_data)

def analyze_contig_sizes(fasta_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    size_data = []
    for root, dirs, files in os.walk(fasta_folder):
        for file in files:
            if file.endswith(".fasta"):
                file_path = os.path.join(root, file)
                sizes = []
                for record in SeqIO.parse(file_path, "fasta"):
                    contig_name = record.id.split('_')[0]
                    size = int(contig_name.split('.')[2][1:])
                    sizes.append(size)
                    size_data.append({"Protein": os.path.splitext(file)[0], "Contig_Size": size})

                # Plot distribution of contig sizes for this file
                plt.figure(figsize=(10, 6))
                sns.histplot(sizes, bins=60, kde=True)  # Augmenter le nombre de bacs (bins) à 60
                plt.title(f"Distribution of Contig Sizes for {os.path.splitext(file)[0]}")
                plt.xlabel("Contig Size")
                plt.ylabel("Frequency")
                plt.tight_layout()
                plot_path = os.path.join(output_folder, f"{os.path.splitext(file)[0]}_contig_size_distribution.png")
                plt.savefig(plot_path)
                plt.close()
    
    # Save combined data to CSV
    df_sizes = pd.DataFrame(size_data)
    df_sizes.to_csv(os.path.join(output_folder, "contig_size_data.csv"), index=False)
    
    # Plot combined distribution of contig sizes
    plt.figure(figsize=(12, 8))
    sns.histplot(df_sizes, x="Contig_Size", hue="Protein", bins=60, kde=True, element="step", common_norm=False)  # Augmenter le nombre de bacs (bins) à 60
    plt.title("Combined Distribution of Contig Sizes")
    plt.xlabel("Contig Size")
    plt.ylabel("Frequency")
    plt.tight_layout()
    combined_plot_path = os.path.join(output_folder, "combined_contig_size_distribution.png")
    plt.savefig(combined_plot_path)
    plt.close()

def main(input_folders, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    all_counts = []
    for folder in input_folders:
        counts = count_sequences(folder)
        counts['Stage'] = os.path.basename(folder)
        all_counts.append(counts)
    
    df = pd.concat(all_counts)
    
    # Save the raw data table
    raw_data_path = os.path.join(output_folder, "orf_counts.csv")
    df.to_csv(raw_data_path, index=False)
    
    # Plot the data
    plt.figure(figsize=(12, 8))
    sns.barplot(x="Protein", y="Count", hue="Stage", data=df)
    plt.xticks(rotation=45)
    plt.title("ORF Counts After Each Filtration Stage")
    plt.tight_layout()

    # Save the plot
    plot_path = os.path.join(output_folder, "orf_counts.png")
    plt.savefig(plot_path)
    plt.close()
    
    # Analyze contig sizes for the third stage only
    analyze_contig_sizes(input_folders[-1], output_folder)
    
    print(f"Analysis completed. Data saved in {output_folder}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python orf_analysis.py <input_folder1> <input_folder2> <input_folder3> <output_folder>")
        sys.exit(1)
    
    input_folders = sys.argv[1:-1]
    output_folder = sys.argv[-1]
    
    main(input_folders, output_folder)
