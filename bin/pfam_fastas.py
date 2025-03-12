
#!/usr/bin/env python    
    
import sys
import os
from Bio import SeqIO
from collections import defaultdict
    
pfam_table = sys.argv[1]
fasta_file = sys.argv[2]
pfam_fams2seqs = defaultdict(set)
    
with open(pfam_table) as fin:
    for line in fin:
        if not line.startswith("#"):
            info = line.strip().split("\t")
            pfam_fams2seqs[info[1]].add(info[0])


seqs2pfam = defaultdict(set)
pfam_small_fams = defaultdict(list)
pfam_singletons = defaultdict()
pfam_fastas = defaultdict(list)

for pfam, seqs in pfam_fams2seqs.items():

    if len(seqs) == 1:
        pfam_singletons[pfam] = seqs

    if len(seqs) == 2:
        pfam_small_fams[pfam] = seqs

    if len(seqs) >= 3:
        pfam_fastas[pfam] = seqs
        
    for s in seqs:
        seqs2pfam[s].add(pfam)


if len(pfam_small_fams) >0:
     with(open('pfam_small_fams.tsv', 'w')) as fout:
        for pfam, seqs in pfam_small_fams.items():
            fout.write('\t'.join([pfam, ','.join(seqs)+'\n']))

if len(pfam_singletons) >0:
    with(open('pfam_singletons.tsv', 'w')) as fout:
        for pfam, seqs in pfam_singletons.items():
            fout.write('\t'.join([pfam, ','.join(seqs)+'\n']))

pfam_clust_mems = open('pfam.clusters_mems.tsv', 'w')
pfam_clust_size = open('pfam.clusters_size.tsv', 'w')
seq2pfam_out = open('pfam_seq2pfam.tsv', 'w')
for pnam, mems in pfam_fams2seqs.items():
    pfam_clust_mems.write('\t'.join([pnam, ','.join(mems)+'\n']))
    pfam_clust_size.write('\t'.join([pnam, str(len(mems))+'\n']))

    for m in mems:
        seq2pfam_out.write('\t'.join([m, pnam+'\n']))

    
pfam_clust_mems.close()
pfam_clust_size.close()




seqs_no_pfam_path = "seqs_no_pfam.faa"
seqs_no_pfam = open(seqs_no_pfam_path, 'w')
    
for record in SeqIO.parse(fasta_file, "fasta"):
    if record.id in seqs2pfam.keys():
        all_pfams = seqs2pfam[record.id]
        for dom in all_pfams:
            
            if dom in pfam_fastas.keys():
                
                path2pfam_fasta =  dom+".pfam.faa"
                with open(path2pfam_fasta, "a") as fout:
                    print(path2pfam_fasta)
                    fout.write(">"+record.id+"\n"+str(record.seq)+"\n")
                    fout.close()
    else:
        seqs_no_pfam.write(">"+record.id+"\n"+str(record.seq)+"\n")
    
seqs_no_pfam.close()