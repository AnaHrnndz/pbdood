#!/usr/bin/env bash
#
# install.sh — Install PBDOOD and all its dependencies.
#
# It does two things:
#   1) Creates (or updates) the conda environment 'pbdood' with all the tools.
#   2) Downloads OG_Delineation pinned to a specific commit into external/
#      (it is not installable via conda/pip, so it is fetched separately).
#
# Usage:
#     bash install.sh
#     conda activate pbdood
#     nextflow run . -profile local --input <fasta> --pfam_db <hmm> --ogd_taxonomy_db <sqlite>
#
set -euo pipefail

# --- Configuration -----------------------------------------------------------
ENV_NAME="pbdood"
OGD_REPO="https://github.com/AnaHrnndz/OG_Delineation.git"
OGD_COMMIT="8ebf3bef818a"          # pinned for reproducibility (2025-11-14)
OGD_DIR="external/OG_Delineation"

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

# --- 1) conda environment ----------------------------------------------------
echo ">> Creating/updating conda environment '${ENV_NAME}'..."
if conda env list | grep -qE "^\s*${ENV_NAME}\s"; then
    "$CONDA" env update -n "$ENV_NAME" -f environment.yml
else
    "$CONDA" env create -f environment.yml
fi

# --- 2) OG_Delineation (external dependency, pinned by commit) ---------------
echo ">> Fetching OG_Delineation @ ${OGD_COMMIT}..."
mkdir -p external
if [ ! -d "${OGD_DIR}/.git" ]; then
    git clone "$OGD_REPO" "$OGD_DIR"
fi
git -C "$OGD_DIR" fetch origin
git -C "$OGD_DIR" checkout "$OGD_COMMIT"

echo
echo ">> Installation complete."
echo "   Activate the environment: conda activate ${ENV_NAME}"
echo "   Then run the pipeline:     nextflow run . -profile local --input <fasta> ..."