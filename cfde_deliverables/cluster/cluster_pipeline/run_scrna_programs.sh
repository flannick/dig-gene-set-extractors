#!/usr/bin/env bash
# scRNA cNMF gene programs -> standard genesets via DIG (geneset-extractors).
#
# Pipeline (4 stages, all DIG-owned):
#   Stage 1: geneset-extractors workflows scrna_cnmf_prepare
#   Stage 2: cnmf prepare / factorize / combine / k_selection_plot  (via generated run_cnmf.sh)
#   Stage 3: cnmf_select_k + cnmf consensus  (via generated run_cnmf_consensus_auto_k.sh)
#   Stage 4: geneset-extractors convert rna_sc_programs  (via generated run_geneset_extractors_from_cnmf.sh)
#
# Usage: run_scrna_programs.sh <matrix.tsv> <meta.tsv> <outdir> <dataset_id>
# Env:
#   CELL_ID_COLUMN      cell ID column in metadata (default: sample_name)
#   CELL_TYPE_COLUMN    cell-type column (default: subclass_label)
#   DONOR_COLUMN        donor column (default: external_donor_name_label)
#   VALUE_TYPE          counts|logcounts (default: counts)
#   ORGANISM            human|mouse (default: human)
#   GENOME_BUILD        hg38|mm10 (default: hg38)
#   K_LIST              K grid for cNMF prepare (default: auto)
#   NHVG                highly-variable genes (default: 2000)
#   NITER               cNMF iterations (default: 100)
#   TOPN                top-N genes per program (default: 100)
#   MAX_CELLS           max cells total (default: 20000)
#   SEED                random seed (default: 1)
#   EXPORT_KIND         score|tpm (default: score)
#   DIG_DIR             path to dig-gene-set-extractors checkout (default: auto-detect)
set -euo pipefail

MATRIX=${1:?matrix_tsv}
META=${2:?meta_tsv}
OUT=${3:?outdir}
NAME=${4:-run}

CID="${CELL_ID_COLUMN:-sample_name}"
CT="${CELL_TYPE_COLUMN:-subclass_label}"
DN="${DONOR_COLUMN:-external_donor_name_label}"
VT="${VALUE_TYPE:-counts}"
ORG="${ORGANISM:-human}"
GBUILD="${GENOME_BUILD:-hg38}"
K_LIST="${K_LIST:-auto}"
NHVG="${NHVG:-2000}"
NITER="${NITER:-100}"
TOPN="${TOPN:-100}"
MAX_CELLS="${MAX_CELLS:-20000}"
SEED="${SEED:-1}"
EXPORT_KIND="${EXPORT_KIND:-score}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Auto-detect DIG dir: either env var, or two levels up from cfde_deliverables/
if [[ -z "${DIG_DIR:-}" ]]; then
  CANDIDATE="$(cd "${HERE}/../../../../" && pwd)"
  if [[ -d "${CANDIDATE}/src/geneset_extractors" ]]; then
    DIG_DIR="${CANDIDATE}"
  else
    echo "DIG_DIR not set and auto-detect failed; set DIG_DIR to dig-gene-set-extractors checkout" >&2
    exit 1
  fi
fi
export PYTHONPATH="${DIG_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}"

PYTHON="${PYTHON_BIN:-python3}"
WORKFLOW_OUT="${OUT}/workflow"
EXTRACTOR_OUT="${OUT}/extractor"
SUBSET_DIR="${WORKFLOW_OUT}/subsets/all"
LOG="${OUT}/run.log"
mkdir -p "${WORKFLOW_OUT}" "${EXTRACTOR_OUT}"

log() { printf '%s\n' "$*" | tee -a "${LOG}"; }

log "[scrna_cnmf 1/4] scrna_cnmf_prepare -> ${WORKFLOW_OUT}"
"${PYTHON}" -m geneset_extractors.cli workflows scrna_cnmf_prepare \
  --matrix_tsv "${MATRIX}" \
  --meta_tsv "${META}" \
  --meta_cell_id_column "${CID}" \
  --cell_type_column "${CT}" \
  --donor_column "${DN}" \
  --matrix_value_type "${VT}" \
  --organism "${ORG}" \
  --genome_build "${GBUILD}" \
  --split_by_cell_type false \
  --max_cells_total "${MAX_CELLS}" \
  --seed "${SEED}" \
  --cnmf_k_list "${K_LIST}" \
  --cnmf_k auto \
  --cnmf_n_iter "${NITER}" \
  --cnmf_numgenes "${NHVG}" \
  --cnmf_export_kind "${EXPORT_KIND}" \
  --cnmf_select_strategy largest_stable \
  --out_dir "${WORKFLOW_OUT}" \
  2>&1 | tee -a "${LOG}"

if [[ ! -d "${SUBSET_DIR}" ]]; then
  log "ERROR: scrna_cnmf_prepare did not produce ${SUBSET_DIR}"; exit 1
fi

log "[scrna_cnmf 2/4] cNMF factorize (cnmf CLI) in ${SUBSET_DIR}"
(cd "${SUBSET_DIR}" && bash run_cnmf.sh 2>&1 | tee -a "${LOG}")

log "[scrna_cnmf 3/4] cNMF consensus (cnmf_select_k + cnmf consensus)"
(cd "${SUBSET_DIR}" && bash run_cnmf_consensus_auto_k.sh 2>&1 | tee -a "${LOG}")

log "[scrna_cnmf 4/4] rna_sc_programs extractor"
(cd "${SUBSET_DIR}" && bash run_geneset_extractors_from_cnmf.sh 2>&1 | tee -a "${LOG}")

# Copy extractor outputs to standard extractor/ location
EXTRACTOR_SRC="$(ls -d "${SUBSET_DIR}/cnmf_out/geneset_extractors_programs_k_"*"_${EXPORT_KIND}" 2>/dev/null | tail -1)"
if [[ -z "${EXTRACTOR_SRC}" ]]; then
  log "ERROR: no rna_sc_programs output found under ${SUBSET_DIR}/cnmf_out/"; exit 1
fi
log "copying ${EXTRACTOR_SRC} -> ${EXTRACTOR_OUT}"
cp -r "${EXTRACTOR_SRC}/." "${EXTRACTOR_OUT}/"

log "[scrna_cnmf] done: ${NAME} -> ${EXTRACTOR_OUT}"
log "  genesets.gmt: $(wc -l < "${EXTRACTOR_OUT}/genesets.gmt" 2>/dev/null || echo '(missing)') programs"
