#!/usr/bin/env bash
#
# install.sh — Install PBDOOD and all its dependencies.
#
# Creates (or updates) the conda environment 'pbdood' with every tool,
# including OG_Delineation, which is now installed via pip from git
# (pinned in environment.yml).
#
# Usage:
#     bash install.sh
#     conda activate pbdood
#     nextflow run . -profile local --input <fasta> --pfam_db <hmm> --ogd_taxonomy_db <sqlite>
#
set -euo pipefail

ENV_NAME="pbdood"

# Use mamba if available (much faster), otherwise fall back to conda.
if command -v mamba >/dev/null 2>&1; then
    CONDA="mamba"
elif command -v conda >/dev/null 2>&1; then
    CONDA="conda"
else
    echo "ERROR: neither 'conda' nor 'mamba' found in PATH." >&2
    echo "Install Miniforge/Miniconda first: https://conda-forge.org/download/" >&2
    exit 1
fi

echo ">> Creating/updating conda environment '${ENV_NAME}'..."
if conda env list | grep -qE "^\s*${ENV_NAME}\s"; then
    "$CONDA" env update -n "$ENV_NAME" -f environment.yml
else
    "$CONDA" env create -f environment.yml
fi

echo
echo ">> Installation complete."
echo "   Activate the environment: conda activate ${ENV_NAME}"
echo "   Then run the pipeline:     nextflow run . -profile local --input <fasta> ..."