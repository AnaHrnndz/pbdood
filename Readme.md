# 🧬 PBDOOD: Phylogeny-based Domain-oriented Orthology Delineation

**PBDOOD** is a bioinformatics pipeline designed to identify **Orthologous Groups (OGs)** across multiple proteomes. It uses the **Nextflow** framework to automate a multi-step process: initial protein clustering based on Pfam domains, phylogenetic tree building, and a final orthology analysis to define OGs.

![Summary of the PBDOOD workflow](images/fig43.png "Summary of the PBDOOD workflow")

---

## ✨ Key Features

* **Domain-based clustering (PfamFams):** Initial clustering groups sequences into Pfam families. Multi-domain sequences are assigned to all relevant Pfam families.
* **De novo clustering (UnkFams):** Sequences without identifiable Pfam domains are clustered separately using **MMseqs2**.
* **Automated gene-tree building:** For families with ≥ 3 sequences, PBDOOD generates multiple sequence alignments (MSAs), removes uninformative columns, and infers a gene tree.
* **Orthology delineation (OGD):** The OGD algorithm analyses each gene tree, detects gene-duplication events, and defines OGs at multiple taxonomic levels, keeping a hierarchical structure.

---

## ⚙️ Pipeline Steps

The pipeline is divided into three main steps:

### 1. Domain and *de novo* clustering
Preparation of the initial protein families (PfamFams and UnkFams).

### 2. Phylogenomics (MSA + trimming + tree building)
Generation of high-quality alignments and inference of gene trees.

### 3. Orthology delineation (OGD)
This algorithm analyses the trees to detect duplications and define OGs. Key steps include:
* **Tree setup:** resolving polytomies, tree rooting, and adding taxonomic annotation from the NCBI taxonomy.
* **Taxonomic misplacement detection:** identification of sequences likely misplaced within the tree.
* **Species-overlap calculation:** for each internal node, the proportion of species shared between the two child nodes relative to the total number of species in the node.
* **Duplication detection:** identification of duplications that give rise to orthologous groups with a low rate of paralogs.

---

## 📋 Requirements and Dependencies

PBDOOD requires:

* **Nextflow** (≥ 23.04) as the workflow manager.
* **conda** or **mamba** to install the rest of the dependencies (recommended), or a container engine (Docker/Apptainer).

All bioinformatics tools are installed automatically through the provided `environment.yml`:

| Category | Tools |
| :--- | :--- |
| **Clustering** | `eggnog-mapper`, `HMMER`, `MMseqs2` |
| **Phylogenomics** | `FAMSA`, `MAFFT`, `FastTree` |
| **Orthology** | `ETE4`, `FastRoot` |

> The orthology step relies on **OG_Delineation** (`og_delineation.py`), which is not distributed via conda/pip. The `install.sh` script fetches it automatically, pinned to a specific commit, into `external/`.

---

## 💻 Installation

### Recommended: conda

1. **Install Nextflow** (requires Java 11+):
   ```bash
   curl -s https://get.nextflow.io | bash
   chmod +x nextflow
   # optionally move it somewhere on your PATH, e.g. ~/.local/bin
   ```
   Alternatively, Nextflow is also available via conda: `conda install -c bioconda nextflow`.

2. **Clone the repository:**
   ```bash
   git clone https://github.com/AnaHrnndz/pbdood.git
   cd pbdood
   ```

3. **Install all dependencies:**
   ```bash
   bash install.sh
   conda activate pbdood
   ```
   This creates the `pbdood` conda environment with every tool and fetches OG_Delineation (pinned) into `external/`.

4. **Download the data bundle:**
   Download the `data` folder from the link below and place it in the repository root (next to `DOOD.nf`):
   `https://saco.csic.es/s/xjzGL86Cj2x2WJs`

   > **Expected project structure:**
   > ```
   > pbdood/
   >     ├── DOOD.nf
   >     ├── nextflow.config
   >     ├── bin/
   >     ├── external/            # OG_Delineation (added by install.sh)
   >     └── data/                # downloaded data bundle
   >         ├── pfam/
   >         ├── ete_taxonomy/
   >         └── Dickeya.fa
   > ```

### Alternative: container

An Apptainer definition file is provided in [`apptainer/`](apptainer/) to build an image with all dependencies. The `docker` and `apptainer` profiles in `nextflow.config` are ready to point at a published image.

---

## 📥 Input Files

* **Proteome FASTA file (required):** a single FASTA with all proteome sequences.
  * **Header convention:** `ncbi_taxid.sequence_name` (e.g. `9606.NP_000001.1`). The species delimiter (`.`) can be changed with `--ogd_sp_delimitator`.
* **Reference databases (required):**
  * **Pfam database** (`Pfam-A.hmm`, v35 by default) — provided in the data bundle under `data/pfam/`.
  * **NCBI taxonomy database** (`e6.taxa.sqlite`) — provided under `data/ete_taxonomy/`.

---

## ▶️ How to Run

Parameters can be set on the command line (`--param value`) or in the config files. Defaults live in `nextflow.config`; `conf/local.config` and `conf/slurm.config` define the executor and resources.

### Local execution

* **Full pipeline:**
  ```bash
  nextflow run . -profile local \
      --input data/Dickeya.fa \
      --pfam_db data/pfam/Pfam-A.hmm \
      --pfam_datadir data \
      --ogd_taxonomy_db data/ete_taxonomy/e6.taxa.sqlite \
      -with-trace -resume
  ```

* **Quick smoke test** (uses the bundled `data/Dickeya.fa`; fill in the DB paths in `conf/test.config` first):
  ```bash
  nextflow run . -profile test,local -resume
  ```

* **Re-run only specific stages** with `-entry`:
  ```bash
  nextflow run . -profile local -entry ogd_rerun -resume
  nextflow run . -profile local -entry summary_only -resume
  ```

* **Run the tree visualisation (smartview):**
  ```bash
  nextflow run . -profile local -entry run_smartview --input_tree path_to_ogd_treefile.nw
  ```

### Cluster execution (SLURM)

Use the `slurm` profile (adjust queue names in `conf/slurm.config` to match your cluster):
```bash
nextflow run . -profile slurm \
    --input <fasta> --pfam_db <hmm> --pfam_datadir <dir> --ogd_taxonomy_db <sqlite> \
    -with-trace -resume
```

---

## 🔧 Key Parameters

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `--input` | *(required)* | Path to the input proteome FASTA. |
| `--pfam_db` | *(required)* | Path to the Pfam HMM database (`Pfam-A.hmm`). |
| `--pfam_datadir` | *(required)* | `data/` folder with Pfam and eggnog-mapper resources. |
| `--ogd_taxonomy_db` | *(required)* | Path to the NCBI taxonomy database (`e6.taxa.sqlite`). |
| `--general_output` | `results/` | Directory where all results are saved. |
| `--hmmer_cpu` | `6` | CPUs for HMMER / hmm_mapper. |
| `--mmseqs_coverage` | `0.3` | Minimum coverage for *de novo* clustering (MMseqs2). |
| `--mmseqs_min_seq_id` | `0.3` | Minimum sequence identity for *de novo* clustering. |
| `--ogd_rooting` | `MinVar` | Tree-rooting algorithm. |
| `--ogd_sp_delimitator` | `.` | Delimiter separating the taxid in FASTA headers. |
| `--ogd_sp_overlap` | `0.1` | Species-overlap threshold for duplication detection. |
| `--ogd_sp_lost` | `0.7` | Species-loss threshold. |

See `nextflow.config` for the full list of parameters.

---

## 📊 Benchmark

We evaluated PBDOOD using the [Quest for Orthologs (QfO)](https://questfororthologs.org/) benchmarking service and the 2022 reference-proteomes dataset, comprising 79 proteomes from Bacteria, Eukaryota and Archaea. All tests are available in the [`benchmark/`](benchmark/) folder; here are some of the results.

![SwissTrees test](benchmark/dood_swisstrees.png "SwissTrees test")
![EC test](benchmark/dood_ec.png "EC test")
![LUCA Generalized Species Tree test](benchmark/dood_std_luca_ncomptrees_rfdist.png "LUCA Generalized Species Tree test")
