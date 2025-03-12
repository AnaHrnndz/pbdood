#!/usr/bin/env nextflow


process pfam_clustering {

    tag { fasta_file }

    publishDir path:  "${params.clustering_output}" , mode:'copy'

    input:
    path fasta_file

    output:
    path "result_hmm_mapper.emapper.hmm_hits", emit: pfam_table

    script:
    """
    python ${params.emapper_dir}/hmm_mapper.py \
        --cut_ga --clean_overlaps clans --usemem \
        --num_servers ${params.hmmer_num_servers} --num_workers ${params.hmmer_num_workers} --cpu ${params.hmmer_cpu} \
        --dbtype hmmdb  -d ${params.emapper_dir}/data/pfam/Pfam-A.hmm \
        --hmm_maxhits 0 --hmm_maxseqlen 60000 \
        --qtype seq -i ${fasta_file} -o result_hmm_mapper  
    """

}

process get_pfam_fastas {

    label 'medium'

    tag { pfam_table }

    memory { params.memory * task.attempt }
    time { params.time * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir path: "${params.fastas_output}" , pattern: "*.pfam.faa", mode: 'copy'
    publishDir path: "${params.clustering_output}", pattern: "pfam_singletons.tsv", mode: 'copy'
    publishDir path: "${params.clustering_output}", pattern: "pfam_small_fams.tsv", mode: 'copy'
    publishDir path: "${params.clustering_output}", pattern: "pfam_seq2pfam.tsv", mode: 'copy'
    publishDir path: "${params.clustering_output}", pattern: "pfam.clusters_size.tsv", mode: 'copy'
    publishDir path: "${params.clustering_output}", pattern: "pfam.clusters_mems.tsv", mode: 'copy'



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
    python ${params.bin_dir}pfam_fastas.py ${pfam_table} ${fasta_file}
    """
    
}

process mmseqs_clustering {

    label 'medium'

    tag { seqs_no_pfam }

    memory { params.memory * task.attempt }
    time { params.time * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir path: "${params.clustering_output}" , mode: 'copy'

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

    python ${params.bin_dir}rename_mmseqs_fams.py ${mmseqs_tsv} ${params.clustering_output}
        
    """
}


process get_mmseqs_fastas {

    label 'medium'

    tag { mmseqs_mems }

    memory { params.memory * task.attempt }
    time { params.time * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2
    
    publishDir path: "${params.fastas_output}", pattern: "mmseqs_*.faa", mode: 'copy'
    publishDir path: "${params.clustering_output}", pattern: "mmseqs_singletons.tsv", mode: 'copy'
    publishDir path: "${params.clustering_output}", pattern: "mmseqs_small_fams.tsv", mode: 'copy'
    publishDir path: "${params.clustering_output}", pattern: "mmseqs_seq2fam.tsv", mode: 'copy'


    input:
    path mmseqs_mems
    path seqs_no_pfam

    output:
    path "mmseqs_*.faa" , emit: mmseqs_fastas
    path "mmseqs_small_fams.tsv", emit: pairs_small_mmseqs
    path "mmseqs_seq2fam.tsv", emit: mmseq_seq2fam
    path "mmseqs_singletons.tsv"


    script:
    """
    python ${params.bin_dir}mmseqs_fastas.py ${mmseqs_mems} ${seqs_no_pfam}
    """
}


workflow CLUSTERING {

    
    fasta_file = Channel.fromPath(params.input)

    
    pfam_clustering(fasta_file)
    get_pfam_fastas(pfam_clustering.out.pfam_table, fasta_file)
    mmseqs_clustering(get_pfam_fastas.out.seqs_no_pfam)
    get_mmseqs_fastas(mmseqs_clustering.out.mmseqs_mems, get_pfam_fastas.out.seqs_no_pfam)


    raw_pfams = get_pfam_fastas.out.pfam_fastas
    raw_mmseqs = get_mmseqs_fastas.out.mmseqs_fastas

    
    all_raw = raw_pfams.merge(raw_mmseqs)

    emit:
    all_raw
}