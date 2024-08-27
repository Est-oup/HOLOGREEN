import os
import sys
import shutil

def create_output_dirs(output_aa_dir, output_aln_dir):
    os.makedirs(output_aa_dir, exist_ok=True)
    os.makedirs(output_aln_dir, exist_ok=True)

def copy_and_rename_files(src_dir, file_pattern, dest_dir, rename_pattern, remove_extension=False, replace_nopolinto=False):
    for file_name in os.listdir(src_dir):
        file_path = os.path.join(src_dir, file_name)
        if os.path.isfile(file_path):
            if file_pattern in file_name:
                base_name = file_name.split(file_pattern)[1]
                if remove_extension:
                    base_name = os.path.splitext(base_name)[0]
                if replace_nopolinto:
                    base_name = base_name.replace("_nopolinto", "_filtered")
                else:
                    base_name = base_name.replace("_polinto", "-polinto_filtered")
                new_name = rename_pattern.format(base_name)
                shutil.copyfile(file_path, os.path.join(dest_dir, new_name))

def main(output_aa_dir, output_aln_dir, alignment_dir, orf_dir, dna_pol_aln_dir, dna_pol_seq_dir, polinto_orf_dir, polinto_aln_dir):
    # Create output directories
    create_output_dirs(output_aa_dir, output_aln_dir)
    
    # Copy and rename ORF sequences
    for file_name in os.listdir(orf_dir):
        if file_name.startswith("db_ORF_pred_derep_") and file_name.endswith("_id_sequences.fasta"):
            base_name = file_name[len("db_ORF_pred_derep_"):-len("_id_sequences.fasta")]
            shutil.copyfile(os.path.join(orf_dir, file_name), os.path.join(output_aa_dir, f"{base_name}.fasta"))

    dna_pol_filtered_seq = os.path.join(dna_pol_seq_dir, "dna_pol_filtered_orfs_sequences.fasta")
    if os.path.isfile(dna_pol_filtered_seq):
        shutil.copyfile(dna_pol_filtered_seq, os.path.join(output_aa_dir, "dnapol_filtered.fasta"))

    for file_name in os.listdir(polinto_orf_dir):
        file_path = os.path.join(polinto_orf_dir, file_name)
        if os.path.isfile(file_path):
            base_name = os.path.splitext(file_name)[0]
            if "_nopolinto" in base_name:
                new_name = base_name.replace("_nopolinto", "_filtered") + ".fasta"
            else:
                new_name = base_name.replace("_polinto", "-polinto_filtered") + ".fasta"
            shutil.copyfile(file_path, os.path.join(output_aa_dir, new_name))

    # Copy and rename alignment files
    copy_and_rename_files(alignment_dir, "db_ORF_pred_derep_", output_aln_dir, "{}", remove_extension=True)
    
    dna_pol_filtered_aln = os.path.join(dna_pol_aln_dir, "dnapol_dna_pol_filtered_orfs_sequences.fasta.m8")
    if os.path.isfile(dna_pol_filtered_aln):
        shutil.copyfile(dna_pol_filtered_aln, os.path.join(output_aln_dir, "dnapol_filtered"))

    for file_name in os.listdir(polinto_aln_dir):
        if file_name.endswith(".log"):
            continue
        file_path = os.path.join(polinto_aln_dir, file_name)
        if os.path.isfile(file_path):
            base_name = os.path.splitext(file_name)[0]
            if "_nopolinto_filtered" in base_name:
                new_name = base_name.replace("_nopolinto_filtered", "_filtered")
            else:
                new_name = base_name.replace("_polinto_filtered", "-polinto_filtered")
            shutil.copyfile(file_path, os.path.join(output_aln_dir, new_name))

if __name__ == "__main__":
    if len(sys.argv) != 9:
        print("Usage: python script.py <output_aa_dir> <output_aln_dir> <alignment_dir> <orf_dir> <dna_pol_aln_dir> <dna_pol_seq_dir> <polinto_orf_dir> <polinto_aln_dir>")
        sys.exit(1)
    output_aa_dir = sys.argv[1]
    output_aln_dir = sys.argv[2]
    alignment_dir = sys.argv[3]
    orf_dir = sys.argv[4]
    dna_pol_aln_dir = sys.argv[5]
    dna_pol_seq_dir = sys.argv[6]
    polinto_orf_dir = sys.argv[7]
    polinto_aln_dir = sys.argv[8]

    main(output_aa_dir, output_aln_dir, alignment_dir, orf_dir, dna_pol_aln_dir, dna_pol_seq_dir, polinto_orf_dir, polinto_aln_dir)
