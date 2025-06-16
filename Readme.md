#AHP 2024


Nextflow pipeline for Clustering, Phylogenomics and Orthology prediction  
  
External tools:  
    clustering module:  
        - eggnog-mapper    
        - HMMER v 3.1b2    
        - BioPython  
        - MMseqs  
    Phylogenomics module:    
        - Famsa 
        - Mafft 
        - Fastree  
    Orthology module:    
        - ETE4  
        - FastRoot (MinVar rooting)
        - Treeprofiler
    
  
nexflow options:  
    -resume  
    -bg > my-file.log  
    -ansi-log false  
    -with-trace 
    -with-timeline qfo_timeline 
    -with-dag flowchart.png

How to run in local:  
    
    - Nexrflow one workflow mode
    bash /data/soft/nextflow run DOOD.nf -c local.config -with-trace -resume
    
    - If you want to repeat just the ogd process, create a new config with the new parametes and the new output folders
    bash /data/soft/nextflow run /data/projects/cpo_pipeline/DOOD.nf -c local.config -with-trace -resume until ogd_pfam,ogd_mmseqs


How to run specific module in local:  
    bash /data/soft/nextflow run subworkflows/phylogenomics.nf -c local.config -resume -entry MODULE_PHYLOGENOMICS  
    bash /data/soft/nextflow run subworkflows/orthology.nf -c local.config -resume -entry MODULE_ORTHOLOGY
    
How to run in cluster:    
    sbatch run_cpo_pipeline.sh cpo_nextflow/DOOD.nf cpo_nextflow/general.config              
  