import pandas as pd

def count_single_contig_clusters(cluster_file):
    """
    Compte le nombre de clusters ayant un seul contig dans un fichier cluster.tsv.
    """
    # Charger le fichier cluster.tsv
    clusters = pd.read_csv(cluster_file, sep="\t", header=None, names=["cluster", "contig"])

    # Affichez les premières lignes pour vérifier
    print("Aperçu des données chargées :")
    print(clusters.head())

    # Vérifier les colonnes uniques
    print(f"Nombre unique de contigs : {clusters['contig'].nunique()}")
    print(f"Nombre unique de clusters : {clusters['cluster'].nunique()}")

    # Compter le nombre de contigs dans chaque cluster
    cluster_counts = clusters.groupby("cluster").size()
    
    # Filtrer les clusters contenant un seul contig
    single_contig_clusters = cluster_counts[cluster_counts == 1]
    
    # Retourner le nombre de ces clusters
    return len(single_contig_clusters)

# Exemple d'utilisation
if __name__ == "__main__":
    cluster_file = "out/bin/clusters.tsv"  # Remplacez par le chemin vers votre fichier cluster.tsv
    single_contig_count = count_single_contig_clusters(cluster_file)
    print(f"Nombre de clusters contenant un seul contig : {single_contig_count}")
