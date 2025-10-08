#!/usr/bin/env nextflow


    
bin_path = "${baseDir}/bin/"
viz_path = "${baseDir}/dood_viz/"
    


process create_output {

    label 'fast'

    publishDir path: { "${params.general_output}" }, mode:'copy'

    output:
    path "clustering/fastas/"
    path "phylogenomics/aln/"
    path "phylogenomics/trim_aln/"
    path "phylogenomics/trees/"
    path "orthology/"
    path "annotation/trees/"
    path "summary"


    script:
    """
    mkdir -p clustering/fastas/
    mkdir -p phylogenomics/aln/
    mkdir -p phylogenomics/trim_aln/
    mkdir -p phylogenomics/trees/
    mkdir -p orthology/
    mkdir -p annotation/trees/
    mkdir -p summary/
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
    publishDir path: "${params.general_output}/clustering/", pattern: "pfam_seq2pfam_info.tsv", mode: 'copy'  
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
    path "pfam_seq2pfam_info.tsv", emit: seq2dom_arq
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
    if [ ${task.attempt} -gt 1 ]
    then    
        FastTree -fastest  ${mmseqs_trim} > ${fasta_name}.nw  
    else
        FastTree ${mmseqs_trim} > ${fasta_name}.nw  
    si
    """
}

process ogd_pfam {

    label 'medium'

    tag { pfam_nw }

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
    path pfam_nw

    output:
    path "*.tree_annot.nw", emit: pfam_ogd_tree
    path "*.ogs_info.tsv", emit: pfam_ogd_info
    path "*.seq2ogs.tsv", emit: pfam_ogd_seq2og 
    path "*.pairs.tsv", emit: pfam_ogd_pairs 
    path "*.stric_pairs.tsv", emit: pfam_ogd_strict_pairs 

    script:
    fasta_name = pfam_nw.baseName
    """
    #mkdir -p ${params.general_output}/orthology/${fasta_name}/

    mkdir -p orthology/${fasta_name}/
    og_delineation.py --tree ${pfam_nw} --output_path ./  \
        --rooting ${params.ogd_rooting} --user_taxonomy ${params.ogd_taxonomy_db} --sp_delimitator  ${params.ogd_sp_delimitator} \
        --sp_ovlap_all ${params.ogd_sp_overlap} --species_losses_perct ${params.ogd_sp_lost} 

    """


}


process ogd_mmseqs {

    label 'medium'

    tag { mmseqs_nw }

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
    path mmseqs_nw

    output:
    path "*.tree_annot.nw", emit: mmseqs_ogd_tree
    path "*.ogs_info.tsv", emit: mmseqs_ogd_info
    path "*.seq2ogs.tsv", emit: mmseqs_ogd_seq2og 
    path "*.pairs.tsv", emit: mmseqs_ogd_pairs 
    path "*.stric_pairs.tsv", emit: mmseqs_ogd_strict_pairs 

    script:
    fasta_name = mmseqs_nw.baseName
    """
    # mkdir -p ${params.general_output}/orthology/${fasta_name}/

    mkdir -p orthology/${fasta_name}/
    og_delineation.py --tree ${mmseqs_nw} --output_path ./  \
        --rooting ${params.ogd_rooting} --user_taxonomy ${params.ogd_taxonomy_db} --sp_delimitator  ${params.ogd_sp_delimitator} \
        --sp_ovlap_all ${params.ogd_sp_overlap} --species_losses_perct ${params.ogd_sp_lost} 

    """


}


process emapper {

    label 'medium'

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

    publishDir path:  "${params.general_output}/annotation/" , mode:'copy'

    input:
    path fasta_file
    

    output:
    path "result_fannot.emapper.annotations", emit: pfam_table_annotation

    script:
    """
    emapper.py --sensmode fast --cpu ${params.hmmer_cpu}  --data_dir ${params.emapper_dir}  \
        -i ${fasta_file} -o result_fannot --output_dir ./
    """

}



process add_func_annot_pfam {

    label 'fast'
    
    tag { pfam_ogd_tree }

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

    publishDir path: { "${params.general_output}/annotation/trees/" }, mode: 'copy'

    input:
    // Recibe una tupla con los 3 archivos
    tuple path(pfam_ogd_tree), path(pfam_table_annotation), path(seq2dom_arq)

    // path pfam_table_annotation
    // path seq2dom_arq
    // path pfam_ogd_tree

    output:
    path "*.fannot.nw", emit: pfam_fannot_tree

    script:
    """
    python ${bin_path}add_annotations.py ${seq2dom_arq} ${pfam_table_annotation} ${pfam_ogd_tree}
    """

}


process add_func_annot_mmseqs {

    label 'fast'
    
    tag { mmseqs_ogd_tree }

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

    publishDir path: { "${params.general_output}/annotation/trees/" }, mode: 'copy'

    input:
    // Recibe una tupla con los 3 archivos
    tuple path(mmseqs_ogd_tree), path(pfam_table_annotation), path(seq2dom_arq)


    output:
    path "*.fannot.nw", emit: mmseqs_fannot_tree

    script:
    """
    python ${bin_path}add_annotations.py ${seq2dom_arq} ${pfam_table_annotation} ${mmseqs_ogd_tree}
    """

}

process run_summary {

    label 'fast'

    publishDir path: { "${params.general_output}/summary/" }, mode: 'copy'

    input:
    val finished_signal
    path fasta_file


    output:
    path 'sp_vs_sp.tsv', emit: sp_vs_sp
    path 'singlecopy_genes.tsv', emit: singlecopy_genes
    path 'ogs_ordered_by_taxid.tsv', emit: ogs_ordered_by_taxid
    path 'dups_counter.tsv', emit: dups_counter
    path 'total_ogs.tsv', emit: total_ogs
    path 'total_pairs.tsv', emit:total_pairs

    script:
    """ 
    python ${bin_path}summary.py ${fasta_file} ${params.general_output}/orthology ${params.ogd_sp_delimitator}
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

    // SUMMARY //

    def pfam_finished = ogd_pfam.out.pfam_ogd_tree.collect()
    def mmseqs_finished = ogd_mmseqs.out.mmseqs_ogd_tree.collect()
    
    
    def sync_channel = pfam_finished.combine(mmseqs_finished)

    // SUMMARY //
    run_summary(sync_channel, fasta_file)

    // ANNOTATION //
    emapper(fasta_file)
    
    ogd_pfam.out.pfam_ogd_tree.combine(emapper.out.pfam_table_annotation).combine(get_pfam_fastas.out.seq2dom_arq).set { pfam_trees_to_annotate_channel }
    add_func_annot_pfam(pfam_trees_to_annotate_channel)

    ogd_mmseqs.out.mmseqs_ogd_tree.combine(emapper.out.pfam_table_annotation).combine(get_pfam_fastas.out.seq2dom_arq).set { mmseqs_trees_to_annotate_channel }   
    add_func_annot_mmseqs(mmseqs_trees_to_annotate_channel)


}

