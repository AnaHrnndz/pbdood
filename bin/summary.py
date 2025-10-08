import glob
from collections import defaultdict
import os
import numpy as np
import pandas as pd
from ete4 import NCBITaxa
import json
import sys

ncbi = NCBITaxa('/data/projects/cpo_pipeline/data/ete_taxonomy/e6.taxa.sqlite')


def get_total_num_seqs_pure_python(fasta_file):
    try:
        count = 0
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
        return count
    except FileNotFoundError:
        print(f"Error: El archivo '{fasta_file}' no se encontró.")
        return 0


def get_ogs_info(list_ogs, output_dir):
        
    
    

    
    total_ogs = list()
    seqs_in_ogs = set()
    ogs_order_by_tax = defaultdict()
    dups_by_tax =  defaultdict(int)
    single_copy = set()
    
    for path2og in list_ogs:
        with open(path2og) as fin:
            for line in fin:
                if line.startswith('#'):
                    continue
                info = line.strip().split('\t')

                ogname = info[0]

                if ogname.endswith('*'):
                    continue

                total_ogs.append(info)

                tax = ogname.split('@')[1].split('|')[0]
                sp_set = set([s.split('.')[0] for s in info[14].split(',')])

                seqs_in_ogs.update(set(info[14].split(',')))

                if len(sp_set) == len(info[14].split(',')):
                    single_copy.add(ogname)

                dups_by_tax[tax] += 1

                lin_set = set()
                for sp in sp_set:
                    lin_set.update(set(ncbi.get_lineage(sp)))


                for l in lin_set:
                    if l not in ogs_order_by_tax.keys():
                        ogs_order_by_tax[l] = defaultdict(list)

                    ogs_order_by_tax[l][tax].append(ogname)

                
    
    
    return total_ogs, seqs_in_ogs, ogs_order_by_tax, dups_by_tax, single_copy
                

def get_pairs_info(list_pairs, sp_delimiter):

    
    sp_vs_sp_counts = defaultdict(lambda: defaultdict(int))
    processed_pairs = set()
    all_species = set()
    total_pairs = set()
    

    for path2pairs in list_pairs:
        name = os.path.basename(path2pairs).replace('.stric_pairs.tsv', '')
        with open(path2pairs) as fin:
            for line in fin:
                info = line.strip().split('\t')
                
                

                sp1 = info[0].split(sp_delimiter)[0]
                sp2 = info[1].split(sp_delimiter)[0]

                if sp1 == sp2:
                    continue

                unique_pair = tuple(sorted((info[0], info[1])))

                # Verificamos si este par ya ha sido procesado.
                if unique_pair not in processed_pairs:
                    
                    all_species.add(sp1)
                    all_species.add(sp2)

                    # Si el par es nuevo, lo agregamos al conjunto y aumentamos el contador.
                    processed_pairs.add(unique_pair)
                    info.append(name)
                    total_pairs.add(tuple(info))
                    sp_vs_sp_counts[sp1][sp2] += 1
                    sp_vs_sp_counts[sp2][sp1] += 1

    species_list = sorted(list(all_species))
    matrix_data = np.zeros((len(species_list), len(species_list)))
    
    for i, sp1 in enumerate(species_list):
        for j, sp2 in enumerate(species_list):
            if sp1 != sp2:
                matrix_data[i, j] = sp_vs_sp_counts[sp1][sp2]

    summary_df = pd.DataFrame(matrix_data, index=species_list, columns=species_list)

    return summary_df, total_pairs



def writint_outputs(output_dir, summary_df, total_pairs, total_ogs, ogs_order_by_tax, dups_by_tax, single_copy):

    summary_outfile = output_dir+'sp_vs_sp.tsv'
    single_copy_outfile = open(output_dir+'singlecopy_genes.tsv', 'w')
    ogs_ordered_outfile = open(output_dir+'ogs_ordered_by_taxid.tsv', 'w')
    dups_counter_outfile = open(output_dir+'dups_counter.tsv', 'w')
    total_ogs_outfile = open(output_dir+'total_ogs.tsv', 'w')
    total_pairs_outfile = open(output_dir+'total_pairs.tsv', 'w') 

    summary_df.to_csv(summary_outfile, sep='\t', index=True, header=True)
    print(f"Tabla de resumen guardada en '{summary_outfile}'")

    single_copy_outfile.write('\n'.join(single_copy))

    for p in total_pairs:
        total_pairs_outfile.write('\t'.join(p)+'\n')

    for tax, ogs_ordered in ogs_order_by_tax.items():
        data = {
            "taxid": tax,
            "ogs": ogs_ordered 
        }
        ogs_ordered_outfile.write(json.dumps(data) + '\n')
    
    for tax, count in dups_by_tax.items():
        dups_counter_outfile.write('\t'.join([str(tax), str(count)+'\n']))

    for og in total_ogs:
        total_ogs_outfile.write('\t'.join(og)+'\n')



    single_copy_outfile.close()
    total_pairs_outfile.close()
    ogs_ordered_outfile.close()
    dups_counter_outfile.close()
    total_ogs_outfile.close()






path2fastafile = sys.argv[1]
total_genes = get_total_num_seqs_pure_python(path2fastafile)

orthology_dir = sys.argv[2]
output_dir = os.getcwd()+'/'

all_ogs = glob.glob(orthology_dir+'/*/'+'*.ogs_info.tsv')

all_seqs2ogs = glob.glob(orthology_dir+'/*/'+'*.seq2ogs.tsv')
all_pairs =  glob.glob(orthology_dir+'/*/'+'*.stric_pairs.tsv')
sp_delimiter = sys.argv[3]


total_ogs, seqs_in_ogs, ogs_order_by_tax, dups_by_tax, single_copy = get_ogs_info(all_ogs, output_dir)

summary_df, total_pairs = get_pairs_info(all_pairs, sp_delimiter)

writint_outputs(output_dir, summary_df, total_pairs, total_ogs, ogs_order_by_tax, dups_by_tax, single_copy)

