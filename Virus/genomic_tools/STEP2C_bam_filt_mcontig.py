import os
import pysam
from multiprocessing import Pool
from tqdm import tqdm


def extract_and_sort_bam(input_bam, temp_bam):
    """
    Extraire les reads multimappés et trier le fichier BAM temporaire.
    """
    temp_bam_sorted = temp_bam.replace(".bam", ".sorted.bam")

    print(f"Extracting and sorting multimapped reads from: {input_bam}")
    command = f"samtools view -@ 21 -b -h -f 0x100 {input_bam} | samtools sort -@ 21 -o {temp_bam_sorted}"
    os.system(command)

    return temp_bam_sorted


def analyze_multimapped_reads(temp_bam):
    """
    Analyser les reads multimappés et déterminer lesquels conserver ou exclure.
    """
    print(f"Analyzing multimapped reads from: {temp_bam}")
    read_dict = {}
    undesirable_reads = set()  # Liste des reads à exclure (flag 0x800)

    with pysam.AlignmentFile(temp_bam, "rb") as bamfile:
        for read in tqdm(bamfile, desc="Analyzing reads"):
            if read.is_unmapped:
                continue

            # Ajouter les reads marqués 0x800 dans la liste des indésirables
            if read.flag & 0x800:
                undesirable_reads.add(read.query_name)
                continue

            # Clé : read_name
            key = read.query_name

            if key not in read_dict:
                read_dict[key] = []
            
            read_dict[key].append(read)

    # Appliquer les règles de filtrage
    filtered_reads = {}
    for key, alignments in read_dict.items():
        if key in undesirable_reads:
            continue  # Ne pas conserver les reads déjà marqués comme indésirables

        contigs = set(aln.reference_name for aln in alignments)
        if len(contigs) > 1:
            # Si le read est mappé sur plusieurs contigs, conserver tout
            filtered_reads[key] = alignments
        else:
            # Si le read est mappé plusieurs fois sur le même contig
            best_alignment = max(alignments, key=lambda x: x.mapping_quality)
            filtered_reads[key] = [best_alignment]

    return undesirable_reads


def filter_chunk(chunk_reads, input_bam, output_chunk, undesirable_reads):
    """
    Filtrer un chunk du fichier BAM en fonction des indésirables.
    """
    with pysam.AlignmentFile(input_bam, "rb") as bamfile, \
         pysam.AlignmentFile(output_chunk, "wb", template=bamfile) as outfile:
        for read in tqdm(chunk_reads, desc=f"Filtering chunk {output_chunk}"):
            if read.query_name not in undesirable_reads:
                outfile.write(read)


def filter_and_write_parallel(input_bam, temp_bam, output_bam, undesirable_reads, num_threads=4):
    """
    Filtrer et écrire le fichier BAM en parallèle.
    """
    print(f"Filtering and writing final BAM for: {input_bam} in parallel")
    with pysam.AlignmentFile(temp_bam, "rb") as bamfile:
        all_reads = list(bamfile)

    # Diviser les reads en chunks
    chunk_size = len(all_reads) // num_threads
    chunks = [all_reads[i:i + chunk_size] for i in range(0, len(all_reads), chunk_size)]

    # Préparer les chemins des chunks de sortie
    output_chunks = [output_bam.replace(".bam", f".chunk{i}.bam") for i in range(len(chunks))]

    # Filtrer les chunks en parallèle
    with Pool(num_threads) as pool:
        pool.starmap(
            filter_chunk,
            [(chunk, temp_bam, output_chunk, undesirable_reads) for chunk, output_chunk in zip(chunks, output_chunks)]
        )

    # Combiner les fichiers BAM filtrés
    print(f"Merging filtered chunks into final BAM: {output_bam}")
    with pysam.AlignmentFile(output_bam, "wb", template=pysam.AlignmentFile(temp_bam, "rb")) as outfile:
        for output_chunk in tqdm(output_chunks, desc="Merging chunks"):
            with pysam.AlignmentFile(output_chunk, "rb") as chunk_file:
                for read in chunk_file:
                    outfile.write(read)
            os.remove(output_chunk)  # Supprimer le chunk temporaire


def process_single_bam(input_bam, temp_dir, output_dir, num_threads=4):
    """
    Processus complet pour un fichier BAM :
    1. Extraction et tri des reads multimappés.
    2. Analyse des alignements.
    3. Écriture du fichier filtré en parallèle.
    """
    # Chemins des fichiers temporaires et de sortie
    temp_bam = os.path.join(temp_dir, os.path.basename(input_bam))
    filtered_bam = os.path.join(output_dir, os.path.basename(input_bam))

    print(f"Processing file: {input_bam}")

    # Étape 1 : Extraction et tri des reads multimappés
    temp_bam_sorted = extract_and_sort_bam(input_bam, temp_bam)

    # Étape 2 : Analyse et identification des indésirables
    undesirable_reads = analyze_multimapped_reads(temp_bam_sorted)

    # Étape 3 : Filtrage et écriture en parallèle
    filter_and_write_parallel(input_bam, temp_bam_sorted, filtered_bam, undesirable_reads, num_threads)

    # Nettoyer les fichiers temporaires
    os.remove(temp_bam_sorted)
    print(f"Processing complete for: {input_bam}")


def process_all_bams(input_dir, temp_dir, output_dir, num_threads=4):
    """
    Traiter tous les fichiers BAM dans un dossier, un par un.
    """
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    bam_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".bam")]

    for bam_file in tqdm(bam_files, desc="Processing all BAM files"):
        process_single_bam(bam_file, temp_dir, output_dir, num_threads)


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 5:
        print("Usage: python script.py <input_dir> <temp_dir> <output_dir> <num_threads>")
        sys.exit(1)

    input_dir = sys.argv[1]
    temp_dir = sys.argv[2]
    output_dir = sys.argv[3]
    num_threads = int(sys.argv[4])

    process_all_bams(input_dir, temp_dir, output_dir, num_threads)
