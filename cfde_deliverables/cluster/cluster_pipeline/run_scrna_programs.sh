#!/usr/bin/env bash
# scRNA -> cNMF gene programs -> contract gene sets (SELF-CONTAINED). Run AFTER setup_env.sh.
# Uses ONLY the installed `cnmf` + `anndata` via scrna_cnmf_programs.py — NO external CLI/toolchain.
# Usage: run_scrna_programs.sh <matrix_cell_by_gene.tsv> <meta.tsv> <outdir> <name>
# Env: CELL_ID_COLUMN CELL_TYPE_COLUMN DONOR_COLUMN MAX_CELLS K NHVG NITER TOPN ORGANISM CITATION
set -uo pipefail
MATRIX=${1:?matrix_tsv}; META=${2:?meta_tsv}; OUT=${3:?outdir}; NAME=${4:-run}
ENV=${ENV:-gsx310}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV" 2>/dev/null || conda activate "$HOME/.conda/envs/$ENV" || { echo "cannot activate $ENV — run setup_env.sh first"; exit 1; }
mkdir -p "$OUT"

echo "[run_scrna] cNMF gene programs for $NAME -> $OUT"
MATRIX_TSV="$MATRIX" META_TSV="$META" OUTDIR="$OUT" NAME="$NAME" \
CELL_ID_COL="${CELL_ID_COLUMN:-cell_id}" CELL_TYPE_COL="${CELL_TYPE_COLUMN:-cell_type}" DONOR_COL="${DONOR_COLUMN:-donor_id}" \
MAX_CELLS="${MAX_CELLS:-20000}" K="${K:-10}" NHVG="${NHVG:-2000}" NITER="${NITER:-20}" \
TOPN="${TOPN:-100}" SEED="${SEED:-14}" ORGANISM="${ORGANISM:-human}" CITATION="${CITATION:-$NAME}" \
  python "$HERE/scrna_cnmf_programs.py" || { echo "cNMF pipeline FAILED — send the error output"; exit 1; }

echo "[run_scrna] package -> SDG submissions (zip + md5)"
[ -d "$OUT/genesets" ] && bash "$HERE/package_submission.sh" "$OUT/genesets" "scrna_cnmf_programs_${NAME}" \
  || echo "  (no $OUT/genesets produced — check the cNMF output above)"
