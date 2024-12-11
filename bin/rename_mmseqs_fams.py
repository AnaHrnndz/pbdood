import random
import string
import sys
from collections import defaultdict


def generate_mmseqs_code():
    # Define the character set: uppercase letters and digits
    chars = string.ascii_uppercase + string.digits
    
    # Generate a random 4-character string
    unique_part = ''.join(random.choices(chars, k=6))
    
    # Combine with the prefix
    mmseqs_code = f"mmseqs_{unique_part}"
    
    return mmseqs_code





mmseqs_tsv = sys.argv[1]
out_dir = sys.argv[2]

mmseqs2seqs = defaultdict(list)
with open(mmseqs_tsv) as fin:
    for line in fin:
        info = line.strip().split('\t')
        mmseqs2seqs[info[0]].append(info[1])
        

mmseqs_mems = open('mmseqs.clusters_mem.tsv', 'w')
mmseqs_size = open('mmseqs.clusters_size.tsv', 'w')
oriname2code = open('mmseqs.ori2code.tsv', 'w')

for original_name, seqs in mmseqs2seqs.items():
    mmseqs_code = generate_mmseqs_code()
    size = str(len(seqs))
    oriname2code.write('\t'.join([original_name, mmseqs_code+'\n']))

    mmseqs_mems.write('\t'.join([mmseqs_code, ','.join(seqs)+'\n']))
    mmseqs_size.write('\t'.join([mmseqs_code, size+'\n']))

mmseqs_mems.close()
mmseqs_size.close()
oriname2code.close()