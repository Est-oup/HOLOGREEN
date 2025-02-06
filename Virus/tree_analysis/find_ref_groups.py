import os
from collections import defaultdict

def read_fasta_ids(fasta_file):
    ids = []
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                ids.append(line[1:].strip())
    return ids

def group_ids_by_prefix(ids):
    grouped_ids = defaultdict(list)
    for id in ids:
        prefix = '_'.join(id.split('_')[:2])
        grouped_ids[prefix].append(id)
    return grouped_ids

def read_fasta_directory(directory):
    all_ids = []
    for root, _, files in os.walk(directory):
        for filename in files:
            file_path = os.path.join(root, filename)
            all_ids.extend(read_fasta_ids(file_path))
    return all_ids

def main(directory):
    all_ids = read_fasta_directory(directory)
    grouped_ids = group_ids_by_prefix(all_ids)
    groups = dict(grouped_ids)
    print("groups = {")
    for prefix, ids in groups.items():
        print(f'    "{prefix}": {ids},')
    print("}")

if __name__ == "__main__":
    directory = 'F:/Virus\\NCLDV\\arbres_rigoux\\seq_references\\'
    main(directory)
