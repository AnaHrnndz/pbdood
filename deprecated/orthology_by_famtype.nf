process ogd_pfam {

    label 'medium'

    tag { pfam_nw }

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
    path pfam_nw

    output:
    path "*.tree_annot.nw", emit: ogd_pfam_tree
    path "*.ogs_info.tsv", emit: ogd_pfam_info
    path "*.seq2ogs.tsv", emit: ogd_pfam_seq2og 
    path "*.pairs.tsv", emit: ogd_pfam_pairs 

    script:
    fasta_name = pfam_nw.baseName
    """
    mkdir -p ${params.orthology_output}/${fasta_name}
    ${params.python_version} ${params.ogd_dir}og_delineation.py --tree ${pfam_nw} --output_path ./  \
        --rooting ${params.ogd_rooting} --user_taxonomy ${params.ogd_taxonomy_db} --sp_delimitator  ${params.ogd_sp_delimitator}
    """


}


process ogd_mmseqs {

    label 'medium'

    tag { mmseqs_nw }

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
    path mmseqs_nw

    output:
    path "*.tree_annot.nw", emit: ogd_mmseqs_tree 
    path "*.ogs_info.tsv", emit: ogd_mmseqs_info
    path "*.seq2ogs.tsv", emit: ogd_mmseqs_seq2og
    path "*.pairs.tsv", emit: ogd_mmseqs_pairs 


    script:
    fasta_name = mmseqs_nw.baseName
    """
    mkdir  -p ${params.orthology_output}/${fasta_name}
    ${params.python_version} ${params.ogd_dir}og_delineation.py --tree ${mmseqs_nw} --output_path ./   \
        --rooting ${params.ogd_rooting} --user_taxonomy ${params.ogd_taxonomy_db} --sp_delimitator  ${params.ogd_sp_delimitator}     
    """

}