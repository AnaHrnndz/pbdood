#!/usr/bin/env nextflow



process align_general {

    label 'medium'

    tag { raw_fasta }

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

    publishDir path: { "${params.phylogenomics_output}/aln/" }, mode: 'copy'

    input:
    path raw_fasta 

    output:
    path "*.aln", emit: general_aln
        
    script:
    fasta_name = raw_fasta.baseName
    """
    num_seqs=\$(grep -c '^>' ${raw_fasta})
    if [ "\$num_seqs" -gt 1000 ]; then
        echo famsa
        famsa -t ${params.famsa_threads} ${raw_fasta} ${fasta_name}.aln 2> align.err  
    else
        echo mafft
        mafft --auto --anysymbol  ${raw_fasta} > ${fasta_name}.aln
    fi
    """

}


process trimming_general {

    label 'fast'

    tag { general_aln }

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

    publishDir path: { "${params.phylogenomics_output}/trim_aln/" }, mode: 'copy'

    input:
    path general_aln

    output:
    path "*.trim", emit: general_trim

    script:
    fasta_name = general_aln.baseName
    """
    ${params.python_version} ${params.bin_dir}trim_alg_v2.py -i ${general_aln} --min_res_abs 3 --min_res_percent 0.1 -o ${fasta_name}.trim
    """

}


process tree_general {

    label 'medium'

    tag { general_trim }

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

    publishDir path: { "${params.phylogenomics_output}/trees/" }, mode: 'copy'

    input:
    path general_trim

    output:
    path "*.nw", emit: general_nw


    script:
    fasta_name = general_trim.baseName
    """
    if [ ${task.attempt} -gt 1 ]
    then    
        FastTree -fastest  ${general_trim}  > ${fasta_name}.nw  
    else 
        FastTree ${general_trim} > ${fasta_name}.nw 
    fi
    """
}









workflow PHYLOGENOMICS {
    
    take:   
    all_raw 
    
    main:
    align_general(all_raw)
    trimming_general(align_general.out.general_aln)
    tree_general(trimming_general.out.general_trim)
       
    emit:
    tree = tree_general.out.general_nw
}

workflow MODULE_PHYLOGENOMICS {
    
    all_raw = Channel.fromPath("${params.fastas_output}*.faa")
                   
    main:
    align_general(all_raw)
    trimming_general(align_general.out.general_aln)
    tree_general(trimming_general.out.general_trim)
    
    emit:
    tree = tree_general.out.general_nw
}