import sys
import os
from ete4 import  PhyloTree
from ete4.smartview import Layout, explorer
from import_layouts import all_layouts



props_popup = set( ["node_name", "dist", "support", "taxid", "rank",
                "common_name", "sci_name", "lca_node_name", "lca_node",
                "evoltype_2", "so_score", "og_name", "dup_node_name",
                "dup_score", "node_create_og", "taxo_outlier",
                "long_branch_outlier", "monophyletic_og", "paraphyletic_og",
                "kegg_ko", "kegg_path", "pref_name", "basal_og", "dom_arq" ])
        

t = PhyloTree(open(sys.argv[1]), parser = 0)
tname = os.path.basename(sys.argv[1])

t.explore(name = tname, layouts = all_layouts, show_leaf_name = False , keep_server=True, port=5000, localhost='0.0.0.0' )
       
        