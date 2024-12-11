#!/usr/bin/env python   

from Bio import SeqIO
from collections import defaultdict
import sys


mmseqs_mems = sys.argv[1]
seqs_no_pfam = sys.argv[2]

seqs2mmseqs = defaultdict()    
mmseqs_sigletons_dict = defaultdict()
mmseqs_small_fams = defaultdict(list)
mmseqs_fams = defaultdict(list)


with open(mmseqs_mems) as fin:
    for line in fin:
        if not line.startswith("#"):
            info = line.strip().split("\t")

            if len(info[1].split(',')) ==1:
                mmseqs_sigletons_dict[info[0]] = info[1]
            
            elif len(info[1].split(',')) ==2:
                mmseqs_small_fams[info[0]] = info[1].split(',')

            elif len(info[1].split(',')) >=3:
                mmseqs_fams[info[0]] = info[1].split(',')

            for s in info[1].split(','):
                    seqs2mmseqs[s] = info[0]



with(open('mmseqs_seq2fam.tsv', 'w')) as fout:
    for seq, mmseq_code in seqs2mmseqs.items():  
        fout.write('\t'.join([seq, mmseq_code+'\n']))

if len(mmseqs_small_fams)>0:
    with(open('mmseqs_small_fams.tsv', 'w')) as fout:
        for mmseq_code, seqs in mmseqs_small_fams.items():  
            fout.write('\t'.join([mmseq_code,','.join(seqs)+'\n']))

if len(mmseqs_sigletons_dict)>0:
    with(open('mmseqs_singletons.tsv', 'w')) as fout:
        for mmseq_code, seq in mmseqs_sigletons_dict.items():  
            fout.write('\t'.join([mmseq_code, seq+'\n']))

for record in SeqIO.parse(seqs_no_pfam, "fasta"):
    mmseq_name_fam = seqs2mmseqs[record.id]
    if mmseq_name_fam in mmseqs_fams.keys():

        path2mmseqs_fasta = mmseq_name_fam + ".faa"
        with open(path2mmseqs_fasta, "a") as fout:
            fout.write(">"+record.id+"\n"+str(record.seq)+"\n")
            fout.close()
