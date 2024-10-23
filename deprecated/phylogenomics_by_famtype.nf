process align_pfam {

    label 'medium'

    tag { raw_pfam_fasta }

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
    path raw_pfam_fasta 

    output:
    path "*.pfam.aln", emit: pfam_aln
        
    script:
    fasta_name = raw_pfam_fasta.baseName
    """
    famsa -t ${params.famsa_threads} ${raw_pfam_fasta} ${fasta_name}.aln 2> align.err  
    """

}

process align_mmseqs {

    label 'medium'

    tag { raw_mmseqs_fasta }

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
    path raw_mmseqs_fasta

    output:
    path "*.mmseqs.aln", emit: mmseqs_aln
    
    
    script:
    fasta_name = raw_mmseqs_fasta.baseName
    """
    famsa -t ${params.famsa_threads} ${raw_mmseqs_fasta} ${fasta_name}.aln 2> align.err
    """

}


process trimming_pfam {

    label 'fast'

    tag { pfam_aln }

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
    path pfam_aln

    output:
    path "*pfam.trim", emit: pfam_trim

    script:
    fasta_name = pfam_aln.baseName
    """
    ${params.python_version} ${params.bin_dir}trim_alg_v2.py -i ${pfam_aln} --min_res_abs 3 --min_res_percent 0.1 -o ${fasta_name}.trim
    """

}

process trimming_mmseqs{

    label 'fast'
    
    tag { mmseqs_aln }

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
    path mmseqs_aln

    output:
    path "*.mmseqs.trim", emit: mmseqs_trim

    script:
    fasta_name = mmseqs_aln.baseName
    """
    ${params.python_version} ${params.bin_dir}trim_alg_v2.py -i ${mmseqs_aln} --min_res_abs 3 --min_res_percent 0.1 -o ${fasta_name}.trim
    """

}

process tree_pfam {

    label 'medium'

    tag { pfam_trim }

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
    path pfam_trim

    output:
    path "*.pfam.nw", emit: pfam_nw


    script:
    fasta_name = pfam_trim.baseName
    """
    if [ ${task.attempt} -gt 1 ]
    then    
        FastTree -fastest  ${pfam_trim}  > ${fasta_name}.nw  
    else 
        FastTree ${pfam_trim} > ${fasta_name}.nw 
    fi
    """
}


process tree_mmseqs {

    label 'medium'

    tag { mmseqs_trim }

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
    path mmseqs_trim

    output:
    path "*.mmseqs.nw", emit: mmseqs_nw


    script:
    fasta_name = mmseqs_trim.baseName
    """
    FastTree ${mmseqs_trim} > ${fasta_name}.nw  
    """
}

