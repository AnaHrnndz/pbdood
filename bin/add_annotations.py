#!/usr/bin/env python    
    
import sys
import os
from Bio import SeqIO
from collections import defaultdict, Counter
from ete4 import PhyloTree
import random


def calculate_best_terms(seq_ids, annotation_dicts, annot_type):
        
        
        term_counter = Counter()
        for seq_id in seq_ids:
            
            annotations_for_seq = annotation_dicts[annot_type].get(seq_id, [])
            
            term_counter.update(annotations_for_seq)

    
        if None in term_counter: # Check if None is actually a key
             del term_counter[None]
        
        if term_counter: # Check if term_counter is not empty
            term, count = term_counter.most_common(1)[0]
            
            percentage = round(((count / len(seq_ids)) * 100), 3)
            best_terms = '|'.join([term, str(percentage)])
       

        else:
            best_terms = None

       
        return best_terms



doms_arq = defaultdict()
seq2info = defaultdict(dict)


domtable = sys.argv[1]
annot_table = sys.argv[2]

for line in open(domtable):
    seqname, dom_arq = line.strip().split('\t')
    doms_arq[seqname] = dom_arq


for line in open(annot_table):
    if line.startswith('#'):
        continue
    info = line.strip().split('\t')
    seq_name = info[0]
    eggnog_ogs = info[4]
    for og in eggnog_ogs.split(','):
        
        level = og.split('|')[0].split('@')[1]
        if level in ['2759', '2', '2157'] :
            basal_og = og.split('|')[0].split('@')[0]
            seq2info['basal_og'][seq_name] = [basal_og]
                
    pref_name = info[8]
    kegg_pathway = info[12]
    kegg_ko = info[11]
    
    seq2info['pref_name'][seq_name] = [pref_name]
    seq2info['kegg_path'][seq_name] = kegg_pathway.split(',')
    seq2info['kegg_ko'][seq_name] = kegg_ko.split(',')



all_props = set( ["node_name", "dist", "support", "taxid", "rank",
                "common_name", "sci_name", "lca_node_name", "lca_node",
                "evoltype_2", "so_score", "og_name", "dup_node_name",
                "dup_score", "node_create_og", "taxo_outlier",
                "long_branch_outlier", "monophyletic_og", "paraphyletic_og",
                "kegg_ko", "kegg_path", "pref_name", "basal_og", "dom_arq" ])

t = PhyloTree(open(sys.argv[3]))
final_tree_name = os.path.basename(sys.argv[3]).replace('.tree_annot.nw', '.fannot.nw')
for n in t.traverse():
        
    if n.is_leaf:
        for annot_type in seq2info.keys():
            lannot = seq2info[annot_type].get(n.name)
            if lannot is not None:
                n.add_prop(annot_type, lannot)

        if n.name in doms_arq.keys():
            domains = doms_arq[n.name]
            n.add_prop('dom_arq', domains)

            
        
    leaves = list(n.leaf_names())
    kko_top_term = calculate_best_terms(leaves, seq2info, "kegg_ko")
    kpath_top_term = calculate_best_terms(leaves, seq2info, "kegg_path")
    pname_top_term = calculate_best_terms(leaves, seq2info, "pref_name")
    basal_og_top_term = calculate_best_terms(leaves, seq2info, "basal_og")

    if kko_top_term: 
        n.add_prop('kegg_ko', kko_top_term)
    if kpath_top_term: 
        n.add_prop('kegg_path', kpath_top_term)
    if pname_top_term: 
        n.add_prop('pref_name', pname_top_term)
    if basal_og_top_term: 
        n.add_prop('basal_og', basal_og_top_term)

    
    random_seq_name = random.choice(list(n.leaf_names()))
    random_node = next(t.search_nodes(name=random_seq_name))
    random_node_domains = random_node.props.get('dom_arq', '')
    n.add_prop('dom_arq', random_node_domains)


t.write(outfile=final_tree_name, props=all_props, format_root_node=True)