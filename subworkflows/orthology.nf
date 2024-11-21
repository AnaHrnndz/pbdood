#!/usr/bin/env nextflow
 
 

process ogd {

    label 'medium'

    tag { general_nw }

    memory { params.memory * task.attempt }
    time { params.time * task.attempt }

    errorStrategy {
    if (task.exitStatus in 137..140) {
        // Exponential backoff strategy: delay increases with each retry
        sleep(Math.pow(2, task.attempt) * 200 as long)
        return 'retry'
    } else {
        return 'ignore'
    }
    }
    maxRetries 3


    publishDir path: { "${params.orthology_output}/${fasta_name}/" }, mode: 'copy'

    input:
    path general_nw

    output:
    path "*.tree_annot.nw", emit: ogd_tree
    path "*.ogs_info.tsv", emit: ogd_info
    path "*.seq2ogs.tsv", emit: ogd_seq2og 
    // path "*.pairs.tsv", emit: ogd_pairs 

    script:
    fasta_name = general_nw.baseName
    """
    mkdir -p ${params.orthology_output}/${fasta_name}
    ${params.python_version} ${params.ogd_dir}og_delineation.py --tree ${general_nw} --output_path ./  \
        --rooting ${params.ogd_rooting} --user_taxonomy ${params.ogd_taxonomy_db} --sp_delimitator  ${params.ogd_sp_delimitator} \
        --sp_ovlap_all 0.1 --species_losses_perct 0.7 --skip_get_pairs 
    """

}







workflow ORTHOLOGY{

    take:
    all_trees
   
    main:
    ogd(all_trees)

    emit:
    ogd.out.ogd_tree
 
}


workflow MODULE_ORTHOLOGY{

    all_trees = Channel.fromPath("${params.phylogenomics_output}trees/*.nw")
    
    ogd(all_trees)

}