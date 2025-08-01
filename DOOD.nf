#!/usr/bin/env nextflow


    
bin_path = "${baseDir}/bin/"
    


process create_output {

    label 'fast'

    publishDir path: { "${params.general_output}" }, mode:'copy'

    output:
    path "clustering/fastas/"
    path "phylogenomics/aln/"
    path "phylogenomics/trim_aln/"
    path "phylogenomics/trees/"
    path "orthology/"


    script:
    """
    mkdir -p clustering/fastas/
    mkdir -p phylogenomics/aln/
    mkdir -p phylogenomics/trim_aln/
    mkdir -p phylogenomics/trees/
    mkdir -p orthology/
    """

}


process pfam_clustering {

    tag { fasta_file }

    publishDir path:  "${params.general_output}/clustering/" , mode:'copy'

    input:
    path fasta_file
    

    output:
    path "result_hmm_mapper.emapper.hmm_hits", emit: pfam_table

    script:
    """
    hmm_mapper.py  --cut_ga --clean_overlaps clans --usemem \
        --num_servers ${params.hmmer_num_servers} --num_workers ${params.hmmer_num_workers} --cpu ${params.hmmer_cpu} \
        --dbtype hmmdb  --data_dir ${params.pfam_datadir} -d ${params.pfam_db} \
        --hmm_maxhits 0 --hmm_maxseqlen 60000 \
        --qtype seq -i ${fasta_file} -o result_hmm_mapper  
    """

}

process get_pfam_fastas {

    label 'medium'

    tag { pfam_table }

    memory { params.memory * task.attempt }
    if (params.time) {
        time { params.time * task.attempt }
    }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir path: "${params.general_output}/clustering/fastas/" , pattern: "*.pfam.faa", mode: 'copy'
    publishDir path: "${params.general_output}/clustering/", pattern: "pfam_singletons.tsv", mode: 'copy'
    publishDir path: "${params.general_output}/clustering/", pattern: "pfam_small_fams.tsv", mode: 'copy'
    publishDir path: "${params.general_output}/clustering/", pattern: "pfam_seq2pfam.tsv", mode: 'copy'
    publishDir path: "${params.general_output}/clustering/", pattern: "pfam.clusters_size.tsv", mode: 'copy'
    publishDir path: "${params.general_output}/clustering/", pattern: "pfam.clusters_mems.tsv", mode: 'copy'



    input:
    path pfam_table
    path fasta_file

    output:
    path "*.pfam.faa", emit: pfam_fastas
    path "seqs_no_pfam.faa", emit: seqs_no_pfam
    path "pfam_small_fams.tsv", emit: pairs_small_pfams
    path "pfam.clusters_mems.tsv"
    path "pfam.clusters_size.tsv"
    path "pfam_seq2pfam.tsv"
    path "pfam_singletons.tsv"
    
    script:
    """
    python ${bin_path}pfam_fastas.py ${pfam_table} ${fasta_file}
    """
    
}

process mmseqs_clustering {

    label 'medium'

    tag { seqs_no_pfam }

    memory { params.memory * task.attempt }
    if (params.time) {
        time { params.time * task.attempt }
    }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir path: "${params.general_output}/clustering/" , mode: 'copy'

    input:
    path seqs_no_pfam


    output:
    path "mmseqs.clusters_mem.tsv", emit: mmseqs_mems
    path "mmseqs.clusters_size.tsv"
    path "mmseqs.clusters.tsv"
    path "mmseqs.ori2code.tsv"


    script:
    mmseqs_db = "seqs_no_pfam.mmseqs_db"
    mmseqs_clustering = "seqs_no_pfam.cluster_db"
    mmseqs_tsv = "mmseqs.clusters.tsv"
    mmseqs_mems = "mmseqs.clusters_mem.tsv"
    mmseqs_size = "mmseqs.clusters_size.tsv"
    """
    mkdir mmseqs_tmp

    mmseqs createdb ${seqs_no_pfam} ${mmseqs_db}
        
    mmseqs cluster ${mmseqs_db} ${mmseqs_clustering} mmseqs_tmp --threads ${params.mmseqs_threads} \
        -c ${params.mmseqs_coverage} --min-seq-id ${params.mmseqs_min_seq_id} -s ${params.sensitivity} \
        --cov-mode ${params.cov_mode} --cluster-mode ${params.cluster_mode}
       
    mmseqs createtsv ${mmseqs_db} ${mmseqs_db} ${mmseqs_clustering} ${mmseqs_tsv}

    python ${bin_path}rename_mmseqs_fams.py ${mmseqs_tsv} ${params.general_output}/clustering/
        
    """
}


process get_mmseqs_fastas {

    label 'medium'

    tag { mmseqs_mems }

    memory { params.memory * task.attempt }
    if (params.time) {
        time { params.time * task.attempt }
    }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2
    
    publishDir path: "${params.general_output}/clustering/fastas/", pattern: "*.mmseqs.faa", mode: 'copy'
    publishDir path: "${params.general_output}/clustering/", pattern: "mmseqs_singletons.tsv", mode: 'copy'
    publishDir path: "${params.general_output}/clustering/", pattern: "mmseqs_small_fams.tsv", mode: 'copy'
    publishDir path: "${params.general_output}/clustering/", pattern: "mmseqs_seq2fam.tsv", mode: 'copy'


    input:
    path mmseqs_mems
    path seqs_no_pfam

    output:
    path "*.mmseqs.faa" , emit: mmseqs_fastas
    path "mmseqs_small_fams.tsv", emit: pairs_small_mmseqs
    path "mmseqs_seq2fam.tsv", emit: mmseq_seq2fam
    path "mmseqs_singletons.tsv"


    script:
    """
    python ${bin_path}mmseqs_fastas.py ${mmseqs_mems} ${seqs_no_pfam}
    """
}



process align_pfam {

    label 'medium'

    tag { raw_pfam_fasta }

    memory { params.memory * task.attempt }
    
    if (params.time) {
        time { params.time * task.attempt }
    }

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

    publishDir path: { "${params.general_output}/phylogenomics/aln/" }, mode: 'copy'

    input:
    path raw_pfam_fasta 

    output:
    path "*.pfam.aln", emit: pfam_aln
        
    script:
    fasta_name = raw_pfam_fasta.baseName
    """
    #famsa -t ${params.famsa_threads} ${raw_pfam_fasta} ${fasta_name}.aln 2> align.err  

    num_seqs=\$(grep -c '^>' ${raw_pfam_fasta})
    if [ "\$num_seqs" -gt 1000 ]; then
        echo famsa
        famsa -t ${params.famsa_threads} ${raw_pfam_fasta} ${fasta_name}.aln 2> align.err  
    else
        echo mafft
        mafft --auto --anysymbol  ${raw_pfam_fasta} > ${fasta_name}.aln
    fi

    """

}

process align_mmseqs {

    label 'medium'

    tag { raw_mmseqs_fasta }

    memory { params.memory * task.attempt }
    if (params.time) {
        time { params.time * task.attempt }
    }

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

    publishDir path: { "${params.general_output}/phylogenomics/aln/" }, mode: 'copy'

    input:
    path raw_mmseqs_fasta

    output:
    path "*.mmseqs.aln", emit: mmseqs_aln
    
    
    script:
    fasta_name = raw_mmseqs_fasta.baseName
    """
    #famsa -t ${params.famsa_threads} ${raw_mmseqs_fasta} ${fasta_name}.aln 2> align.err

    num_seqs=\$(grep -c '^>' ${raw_mmseqs_fasta})
    if [ "\$num_seqs" -gt 1000 ]; then
        echo famsa
        famsa -t ${params.famsa_threads} ${raw_mmseqs_fasta} ${fasta_name}.aln 2> align.err  
    else
        echo mafft
        mafft --auto --anysymbol  ${raw_mmseqs_fasta} > ${fasta_name}.aln
    fi
    """

}


process trimming_pfam {

    label 'fast'

    tag { pfam_aln }

    memory { params.memory * task.attempt }
    if (params.time) {
        time { params.time * task.attempt }
    }

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

    publishDir path: { "${params.general_output}/phylogenomics/trim_aln/" }, mode: 'copy'

    input:
    path pfam_aln

    output:
    path "*pfam.trim", emit: pfam_trim

    script:
    fasta_name = pfam_aln.baseName
    """
    python ${bin_path}trim_alg_v2.py -i ${pfam_aln} --min_res_abs 3 --min_res_percent 0.1 -o ${fasta_name}.trim
    """

}

process trimming_mmseqs{

    label 'fast'
    
    tag { mmseqs_aln }

    memory { params.memory * task.attempt }
    if (params.time) {
        time { params.time * task.attempt }
    }

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

    publishDir path: { "${params.general_output}/phylogenomics/trim_aln/" }, mode: 'copy'

    input:
    path mmseqs_aln

    output:
    path "*.mmseqs.trim", emit: mmseqs_trim

    script:
    fasta_name = mmseqs_aln.baseName
    """
    python ${bin_path}trim_alg_v2.py -i ${mmseqs_aln} --min_res_abs 3 --min_res_percent 0.1 -o ${fasta_name}.trim
    """

}

process tree_pfam {

    label 'medium'

    tag { pfam_trim }

    memory { params.memory * task.attempt }
    if (params.time) {
        time { params.time * task.attempt }
    }

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

    publishDir path: { "${params.general_output}/phylogenomics/trees/" }, mode: 'copy'

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
    if (params.time) {
        time { params.time * task.attempt }
    }

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

    publishDir path: { "${params.general_output}/phylogenomics/trees/" }, mode: 'copy'

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

process ogd_pfam {

    label 'medium'

    tag { general_nw }

    memory { params.memory * task.attempt }
    if (params.time) {
        time { params.time * task.attempt }
    }
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


    publishDir path: { "${params.general_output}/orthology/${fasta_name}/" }, mode: 'copy'

    input:
    path general_nw

    output:
    path "*.tree_annot.nw", emit: ogd_tree
    path "*.ogs_info.tsv", emit: ogd_info
    path "*.seq2ogs.tsv", emit: ogd_seq2og 
    path "*.pairs.tsv", emit: ogd_pairs 
    path "*.stric_pairs.tsv", emit: ogd_strict_pairs 

    script:
    fasta_name = general_nw.baseName
    """
    #mkdir -p ${params.general_output}/orthology/${fasta_name}/

    mkdir -p orthology/${fasta_name}/
    og_delineation.py --tree ${general_nw} --output_path ./  \
        --rooting ${params.ogd_rooting} --user_taxonomy ${params.ogd_taxonomy_db} --sp_delimitator  ${params.ogd_sp_delimitator} \
        --sp_ovlap_all ${params.ogd_sp_overlap} --species_losses_perct ${params.ogd_sp_lost} 

    """


}


process ogd_mmseqs {

    label 'medium'

    tag { general_nw }

    memory { params.memory * task.attempt }
    if (params.time) {
        time { params.time * task.attempt }
    }

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


    publishDir path: { "${params.general_output}/orthology/${fasta_name}/" }, mode: 'copy'

    input:
    path general_nw

    output:
    path "*.tree_annot.nw", emit: ogd_tree
    path "*.ogs_info.tsv", emit: ogd_info
    path "*.seq2ogs.tsv", emit: ogd_seq2og 
    path "*.pairs.tsv", emit: ogd_pairs 
    path "*.stric_pairs.tsv", emit: ogd_strict_pairs 

    script:
    fasta_name = general_nw.baseName
    """
    # mkdir -p ${params.general_output}/orthology/${fasta_name}/

    mkdir -p orthology/${fasta_name}/
    og_delineation.py --tree ${general_nw} --output_path ./  \
        --rooting ${params.ogd_rooting} --user_taxonomy ${params.ogd_taxonomy_db} --sp_delimitator  ${params.ogd_sp_delimitator} \
        --sp_ovlap_all ${params.ogd_sp_overlap} --species_losses_perct ${params.ogd_sp_lost} 

    """


}

workflow {

    create_output()
    fasta_file = Channel.fromPath(params.input)
    

    // CLUSTERING //
    pfam_clustering(fasta_file)
    get_pfam_fastas(pfam_clustering.out.pfam_table, fasta_file)
    mmseqs_clustering(get_pfam_fastas.out.seqs_no_pfam)
    get_mmseqs_fastas(mmseqs_clustering.out.mmseqs_mems, get_pfam_fastas.out.seqs_no_pfam)

    raw_pfams = get_pfam_fastas.out.pfam_fastas.flatten()
    raw_mmseqs = get_mmseqs_fastas.out.mmseqs_fastas.flatten()

    // PHYLOGENOMICS //
    align_pfam(raw_pfams)
    trimming_pfam(align_pfam.out.pfam_aln)
    tree_pfam(trimming_pfam.out.pfam_trim)

    align_mmseqs(raw_mmseqs)
    trimming_mmseqs(align_mmseqs.out.mmseqs_aln)
    tree_mmseqs(trimming_mmseqs.out.mmseqs_trim)
    

    // ORTHOLOGY //
    ogd_pfam(tree_pfam.out.pfam_nw)
    ogd_mmseqs(tree_mmseqs.out.mmseqs_nw)  
}