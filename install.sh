#!/usr/bin/env bash
#
# install.sh — Instalación de PBDOOD y todas sus dependencias.
#
# Hace dos cosas:
#   1) Crea (o actualiza) el entorno conda 'pbdood' con todas las herramientas.
#   2) Descarga OG_Delineation fijado a un commit concreto en external/
#      (no es instalable por conda/pip, por eso se obtiene aparte).
#
# Uso:
#     bash install.sh
#     conda activate pbdood
#     nextflow run . -profile local --input <fasta> --pfam_db <hmm> --ogd_taxonomy_db <sqlite>
#
set -euo pipefail

# --- Configuración -----------------------------------------------------------
ENV_NAME="pbdood"
OGD_REPO="https://github.com/AnaHrnndz/OG_Delineation.git"
OGD_COMMIT="8ebf3bef818a"          # fijado para reproducibilidad (2025-11-14)
OGD_DIR="external/OG_Delineation"

# Detecta mamba si está disponible (mucho más rápido), si no usa conda.
if command -v mamba >/dev/null 2>&1; then
    CONDA="mamba"
elif command -v conda >/dev/null 2>&1; then
    CONDA="conda"
else
    echo "ERROR: no se encuentra 'conda' ni 'mamba' en el PATH." >&2
    echo "Instala Miniforge/Miniconda primero: https://conda-forge.org/download/" >&2
    exit 1
fi

# --- 1) Entorno conda --------------------------------------------------------
echo ">> Creando/actualizando el entorno conda '${ENV_NAME}'..."
if conda env list | grep -qE "^\s*${ENV_NAME}\s"; then
    "$CONDA" env update -n "$ENV_NAME" -f environment.yml
else
    "$CONDA" env create -f environment.yml
fi

# --- 2) OG_Delineation (dependencia externa, fijada por commit) --------------
echo ">> Obteniendo OG_Delineation @ ${OGD_COMMIT}..."
mkdir -p external
if [ ! -d "${OGD_DIR}/.git" ]; then
    git clone "$OGD_REPO" "$OGD_DIR"
fi
git -C "$OGD_DIR" fetch origin
git -C "$OGD_DIR" checkout "$OGD_COMMIT"

echo
echo ">> Instalación completada."
echo "   Activa el entorno con:   conda activate ${ENV_NAME}"
echo "   Luego lanza el pipeline: nextflow run . -profile local --input <fasta> ..."