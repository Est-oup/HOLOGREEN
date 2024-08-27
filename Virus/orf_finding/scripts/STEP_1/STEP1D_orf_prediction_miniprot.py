import os
import sys
import subprocess
from Bio import SeqIO

# Vérification des arguments de la ligne de commande
if len(sys.argv) != 6:
    print("Usage: python script.py <info_folder> <query_contigs_folder> <seq_references_folder> <output_folder> <log_folder>")
    sys.exit(1)

info_folder = sys.argv[1]
query_contigs_folder = sys.argv[2]
seq_references_folder = sys.argv[3]
output_folder = sys.argv[4]
log_folder = sys.argv[5]
orfpredminiprot_folder = output_folder  # Pas de sous-dossier supplémentaire

def list_files(info_folder):
    files = os.listdir(info_folder)
    return files

def load_sequences(query_contigs_folder, seq_references_folder):
    query_contigs = {}
    seq_references = {}

    for file in os.listdir(query_contigs_folder):
        file_path = os.path.join(query_contigs_folder, file)
        with open(file_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                query_contigs[record.id.split()[0]] = str(record.seq)

    for file in os.listdir(seq_references_folder):
        protein_name = file.split('.')[0]
        file_path = os.path.join(seq_references_folder, file)
        seq_references[protein_name] = {}
        with open(file_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq_references[protein_name][record.id] = (str(record.seq), file)
    
    return query_contigs, seq_references

def log_message(message):
    with open(os.path.join(log_folder, "log_miniprot.txt"), "a") as log_file:
        log_file.write(message + "\n")
    print(message)

def run_miniprot(contig_seq, protein_seq, output_file):
    query_file = "temp_query.fasta"
    ref_file = "temp_ref.fasta"

    with open(query_file, "w") as qf:
        qf.write(f">query\n{contig_seq}\n")

    with open(ref_file, "w") as rf:
        rf.write(f">ref\n{protein_seq}\n")

    command = f'miniprot --trans {query_file} {ref_file} > "{output_file}"'
    try:
        subprocess.run(command, shell=True, check=True)
        log_message(f"Ran Miniprot on {query_file} against {ref_file}, output: {output_file}")
        if os.path.getsize(output_file) == 0:
            log_message(f"Miniprot output file is empty: {output_file}")
            return None
    except subprocess.CalledProcessError as e:
        log_message(f"Error running Miniprot on {query_file} against {ref_file}: {e}")
        return None
    finally:
        os.remove(query_file)
        os.remove(ref_file)
    return output_file

def get_protein_type_from_filename(filename):
    parts = filename.split('_')
    if len(parts) >= 2:
        return parts[1]
    return "unknown"

def process_file(file, query_contigs, seq_references):
    protein_type = file.split('_')[1]
    file_path = os.path.join(info_folder, file)
    results = []

    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            contig_id = parts[0]
            protein_id = parts[1]
            
            contig_seq = query_contigs.get(contig_id)
            protein_data = seq_references.get(protein_type, {}).get(protein_id)
            if protein_data:
                protein_seq, _ = protein_data
            else:
                protein_seq = None

            if contig_seq and protein_seq:
                output_filename = f"{contig_id}_{protein_id}.gff"
                protein_output_folder = os.path.join(orfpredminiprot_folder, protein_type)
                os.makedirs(protein_output_folder, exist_ok=True)
                output_file = os.path.join(protein_output_folder, output_filename)
                result = run_miniprot(contig_seq, protein_seq, output_file)
                if result:
                    results.append(result)
                else:
                    log_message(f"Miniprot did not produce output for contig: {contig_id}, protein: {protein_id}")
            else:
                if not contig_seq:
                    log_message(f"Séquence non trouvée pour le contig : {contig_id}")
                if not protein_seq:
                    log_message(f"Séquence non trouvée pour la protéine : {protein_id}")
    return results

def main():
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(log_folder, exist_ok=True)

    query_contigs, seq_references = load_sequences(query_contigs_folder, seq_references_folder)
    files = list_files(info_folder)

    total_files = 0
    total_orfs = 0
    errors = 0
    empty_files = 0

    for file in files:
        try:
            results = process_file(file, query_contigs, seq_references)
            total_files += 1
            total_orfs += len(results)
            if not results:
                empty_files += 1
        except Exception as exc:
            log_message(f"{file} a généré une exception: {exc}")
            errors += 1

    # Write log file
    log_file_path = os.path.join(log_folder, "log_miniprot_summary.txt")
    with open(log_file_path, 'w') as log_file:
        log_file.write(f"Total files processed: {total_files}\n")
        log_file.write(f"Total ORFs found: {total_orfs}\n")
        log_file.write(f"Total errors: {errors}\n")
        log_file.write(f"Total empty Miniprot output files: {empty_files}\n")

if __name__ == "__main__":
    main()
