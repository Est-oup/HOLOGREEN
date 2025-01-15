import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def load_paf(file_path):
    """Charge le fichier PAF et valide les colonnes."""
    try:
        paf = pd.read_csv(file_path, sep="\t", header=None, usecols=range(12))
        paf.columns = [
            "query", "query_length", "query_start", "query_end",
            "strand", "target", "target_length", "target_start",
            "target_end", "residue_matches", "alignment_block_length", "mapping_quality"
        ]
        # Convertir les colonnes numériques
        numeric_cols = ["query_length", "query_start", "query_end",
                        "target_length", "target_start", "target_end",
                        "residue_matches", "alignment_block_length", "mapping_quality"]
        paf[numeric_cols] = paf[numeric_cols].apply(pd.to_numeric, errors="coerce")
        
        # Vérifier les valeurs NaN
        if paf[numeric_cols].isnull().any().any():
            logging.warning("Certaines colonnes contiennent des valeurs non numériques.")
            logging.warning("Lignes problématiques :")
            logging.warning(paf[paf[numeric_cols].isnull().any(axis=1)])
        return paf
    except Exception as e:
        logging.error(f"Erreur lors du chargement du fichier PAF : {e}")
        raise

def filter_alignments(paf, min_overlap):
    """Filtre les alignements avec un chevauchement significatif."""
    return paf[
        (paf["query"] != paf["target"]) &  # Pas d'alignement d'une séquence sur elle-même
        (paf["query_end"] - paf["query_start"] > min_overlap) &
        (paf["target_end"] - paf["target_start"] > min_overlap)
    ]

def main():
    paf_file = "minimap.paf"
    min_overlap = 100  # Par exemple, chevauchement minimum de 100 pb

    logging.info("Début du processus.")
    paf = load_paf(paf_file)

    if paf.empty:
        logging.warning("Aucun alignement trouvé dans le fichier PAF.")
        return

    logging.info("Filtrage des alignements significatifs.")
    significant_alignments = filter_alignments(paf, min_overlap)
    logging.info(f"{len(significant_alignments)} alignements significatifs trouvés.")

    if not significant_alignments.empty:
        significant_alignments.to_csv("significant_overlaps.csv", index=False)
        logging.info("Alignements significatifs enregistrés dans 'significant_overlaps.csv'.")
    else:
        logging.info("Aucun alignement significatif trouvé.")

if __name__ == "__main__":
    main()
