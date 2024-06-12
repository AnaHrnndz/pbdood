#AHP 2024


Nextflow pipeline for Clustering, Phylogenomics and Orthology prediction  
  
External tools:  
- eggnog-mapper  
- HMMER v 3.1b2  
- BioPython  
- MMseqs  
- Famsa  
- Fastree  
- ETE4  
    
  
nexflow options:  
    -resume  
    -bg > my-file.log  
    -ansi-log false  

How to run in local:  
    /data/soft/nextflow -C local.config run cpo.nf   
    
How to run in cluster:    
    sbatch run_cpo_pipeline.sh cpo_nextflow/cpo.nf cpo_nextflow/general.config              
  