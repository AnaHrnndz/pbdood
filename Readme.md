
# DOOD: Domain Oriented Orthology Delineation

---

## Overview

**DOOD (Domain Oriented Orthology Delineation)** is a bioinformatics tool that identifies orthologous groups across different proteomes. It first clusters proteins based on Pfam domains, then builds phylogenetic trees, and finally uses the OGD Python script to analyze these trees, pinpointing duplication events and defining the orthologous groups.
## Key Features

* **Domain-Based Clustering:** Initial protein clustering based on Pfam domains, grouping proteins into Pfam families (PfamFams).
* **Handling Multi-Domain:** Multi-domain sequences will be part of every PfamFam they map to
* **Unknown Sequences:** Sequences lacking Pfam domains are clustered using MMseqs to form UnkFams.
* **Automated Gene Tree Building:** For families with at least three sequences, DOOD automatically generates multiple sequence alignments (MSAs), with an in-house script removes uninformative columns from MSAs and finllany infers a gene trees.
* **Orthology Delineation (OGD):** Employs the OGD algorithm to analyze each gene trees, detect duplication events, and define orthologous groups (OGs) at various taxonomic levels.



![Alt text for the image](images/graphical_summary.jpg "Summary")

---

## Requirements

DOOD uses several external bioinformatics tools and libraries. Make sure these are accessible in your environment.

### Workflow Manager

* `Nextflow`

### Clustering Module Dependencies

* `eggnog-mapper`
* `HMMER v 3.1b2`
* `BioPython`
* `MMseqs`

### Phylogenomics Module Dependencies

* `Famsa`
* `Mafft`
* `Fasttree`

### Orthology Module Dependencies

* `ETE4`
* `FastRoot` 
* `Treeprofiler`

---

## Input

DOOD requires the following input files:

* **Proteome FASTA file:** A single FASTA file containing all proteome sequences.
    * **Sequence Naming Convention:** Sequence headers must follow the format `ncbi_taxid.sequence_name`. The `.` can be replaced by other characters like `@` or `|`.
    * **Example:** `9606.NP_000001.1`

### Optional Input Files

* **Pfam database:** If you wish to use a specific version of the Pfam database, provide its path.
* **Species tree:** Provide a custom species tree if required for your analysis.
* **NCBI taxonomy database:** Provide a specific version of the NCBI taxonomy database if needed.

---

## Output

DOOD create 3 main output folder:

* **Clustering:** Here you can find all the output for the clustering steps.
    * fasta: This folder contain all the raw fasta created for each family
    * result_hmm_mapper.emapper.hmm_hits: is the main output file from the pfam clustering
    

* **Phylogenomics:**  
    * aln: folder with all the alignments
    * trim_aln: folder with the trimmed aligments
    * trees: folder with all the trees

* **Orthology:** This folder contains a subfolder per family. The main outputs from OGD are:
    * *.tree_annot.nw
    * *.ogs_info.tsv
---

## How to Run DOOD

DOOD can be executed locally or on an HPC cluster using Nextflow.

### Local Execution

To run DOOD on your local machine:

* **Run the entire pipeline:**
    ```bash
    bash /path/to/nextflow run DOOD.nf -c local.config -with-trace -resume
    ```
    * Replace `/path/to/nextflow` with the actual path to your Nextflow executable.

* **Repeat specific OGD processes:**
    If you need to re-run only the OGD (Orthologs Group Delineation) part, create a new configuration file (`new_local.config`) with your desired parameters and new output directories, then execute:
    ```bash
    bash /path/to/nextflow run /data/projects/cpo_pipeline/DOOD.nf -c new_local.config -with-trace -resume until ogd_pfam,ogd_mmseqs
    ```

* **Run specific modules:**
    For targeted execution of individual modules:
    ```bash
    bash /path/to/nextflow run subworkflows/phylogenomics.nf -c local.config -resume -entry MODULE_PHYLOGENOMICS
    bash /path/to/nextflow run subworkflows/orthology.nf -c local.config -resume -entry MODULE_ORTHOLOGY
    ```

### HPC Execution

To submit DOOD to an HPC cluster using Slurm (or a similar job scheduler):

```bash
sbatch run_cpo_pipeline.sh cpo_nextflow/DOOD.nf cpo_nextflow/general.config