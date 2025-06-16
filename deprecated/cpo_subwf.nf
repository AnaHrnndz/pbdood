#!/usr/bin/env nextflow

/*
    IMPORT SUBWORKFLOWS
*/

include {CLUSTERING     } from './subworkflows/clustering'
include {PHYLOGENOMICS  } from './subworkflows/phylogenomics'
include {ORTHOLOGY  } from './subworkflows/orthology'

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


workflow{

    create_output()

    all_fastas = CLUSTERING()
    
    all_trees = PHYLOGENOMICS(all_fastas.flatten())
    
    ORTHOLOGY(all_trees.flatten())

}