#!/bin/bash
# run_catlas_fetal.sh — fetal CATlas gene sets + adult-fetal maturation contrast pipeline.
# Sources: GSE149683 (fetal; BICCN/NIH; sci-ATAC-seq3; public) + Batch 13 adult sets.
#
# Usage:
#   bash run_catlas_fetal.sh setup    # download File_S6 + File_S2 (fetal metadata)
#   bash run_catlas_fetal.sh peek     # print first 3 lines of File_S6 to inspect format
#   bash run_catlas_fetal.sh run      # derive fetal raw + controlled gene sets
#   bash run_catlas_fetal.sh contrast # derive 6-way maturation contrast sets
#   bash run_catlas_fetal.sh combine  # package -> zip -> copy to SDG
#   bash run_catlas_fetal.sh all      # setup -> run -> contrast -> combine
#
# Env knobs: WORK, ADULT_SETS, SDG, ACT_THRESH (default 0), LOW (default 0.25), MIN_GENES (default 10)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
export WORK="${WORK:-/humgen/diabetes2/users/gage/CFDE/catlas_work}"
FETAL_WORK="${WORK}/fetal"
mkdir -p "$FETAL_WORK"
export FILE_S6="${FETAL_WORK}/GSE149683_File_S6.csv.gz"
FILE_S2="${FETAL_WORK}/GSE149683_File_S2.metadata.txt.gz"
export OUTDIR_FETAL="${FETAL_WORK}/catlas_fetal_genesets_full"
export OUTDIR_CONTRAST="${FETAL_WORK}/catlas_maturation_contrast"
export ADULT_SETS="${ADULT_SETS:-${WORK}/catlas_genesets_full}"
SDG="${SDG:-/humgen/diabetes2/users/ryank/CFDE/geneset_extractors/submissions/SDG}"
GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149683/suppl"

dl(){
  local url="$1" dest="$2"
  [ -s "$dest" ] && { echo "have $(basename "$dest")"; return; }
  echo "downloading $(basename "$dest") ..."
  wget -q -O "$dest.part" "$url" && mv "$dest.part" "$dest" || { echo "FAILED: $url"; exit 1; }
}

setup(){
  dl "${GEO_BASE}/GSE149683_File_S6.Cicero_gene_activity_scores_by_cell_type.csv.gz" "$FILE_S6"
  dl "${GEO_BASE}/GSE149683_File_S2.Metadata_of_high_quality_cells.txt.gz" "$FILE_S2"
  echo "File_S6 size: $(ls -lh "$FILE_S6" | awk '{print $5}')"
  echo "File_S2 size: $(ls -lh "$FILE_S2" | awk '{print $5}')"
}

peek(){
  echo "=== First 3 lines of File_S6 ==="
  zcat "$FILE_S6" | head -3 | cut -c1-200 || true
  echo ""
  echo "=== Column count (line 1) ==="
  zcat "$FILE_S6" | head -1 | tr ',' '\n' | wc -l || true
  echo "=== Row count ==="
  zcat "$FILE_S6" | wc -l || true
}

run_fetal(){
  echo "Deriving fetal gene sets from File_S6 ..."
  mkdir -p "$OUTDIR_FETAL"
  WORK="$WORK" FILE_S6="$FILE_S6" OUTDIR_FETAL="$OUTDIR_FETAL" \
  ACT_THRESH="${ACT_THRESH:-0}" LOW="${LOW:-0.25}" MIN_GENES="${MIN_GENES:-10}" \
    python3 "$HERE/catlas_fetal_genesets.py"
  echo "Fetal sets written: $(ls -d "$OUTDIR_FETAL"/CATLAS_FETAL_*_accessible_raw 2>/dev/null | wc -l) raw, $(ls -d "$OUTDIR_FETAL"/CATLAS_FETAL_*_specific_Up 2>/dev/null | wc -l) controlled"
}

contrast(){
  echo "Deriving maturation contrast sets ..."
  [ -d "$ADULT_SETS" ] || { echo "ERROR: ADULT_SETS not found: $ADULT_SETS"; exit 1; }
  [ -d "$OUTDIR_FETAL" ] || { echo "ERROR: run 'run' phase first"; exit 1; }
  mkdir -p "$OUTDIR_CONTRAST"
  WORK="$WORK" ADULT_SETS="$ADULT_SETS" FETAL_SETS="$OUTDIR_FETAL" \
  OUTDIR_CONTRAST="$OUTDIR_CONTRAST" MIN_GENES="${MIN_GENES:-10}" \
    python3 "$HERE/catlas_maturation_contrast.py"
  echo "Contrast sets: $(ls -d "$OUTDIR_CONTRAST"/CATLAS_* 2>/dev/null | wc -l)"
}

combine(){
  echo "Packaging fetal + contrast gene sets ..."
  local zip_fetal="${FETAL_WORK}/batch_catlas_fetal_genesets.zip"
  local zip_contrast="${FETAL_WORK}/batch_catlas_maturation_contrast.zip"
  local n_fetal n_contrast
  n_fetal=$(ls -d "$OUTDIR_FETAL"/CATLAS_FETAL_* 2>/dev/null | wc -l)
  n_contrast=$(ls -d "$OUTDIR_CONTRAST"/CATLAS_* 2>/dev/null | wc -l)
  echo "fetal sets: $n_fetal | contrast sets: $n_contrast"
  ( cd "$OUTDIR_FETAL/.." && zip -qr "$zip_fetal" "$(basename "$OUTDIR_FETAL")" )
  ( cd "$OUTDIR_CONTRAST/.." && zip -qr "$zip_contrast" "$(basename "$OUTDIR_CONTRAST")" )
  if [ -d "$SDG" ]; then
    cp "$zip_fetal" "$SDG/"; echo "copied fetal -> $SDG"
    cp "$zip_contrast" "$SDG/"; echo "copied contrast -> $SDG"
  fi
  md5sum "$zip_fetal" "$zip_contrast"
}

case "${1:-all}" in
  setup)    setup;;
  peek)     peek;;
  run)      run_fetal;;
  contrast) contrast;;
  combine)  combine;;
  all)      setup; peek; run_fetal; contrast; combine;;
  *) echo "usage: $0 {setup|peek|run|contrast|combine|all}"; exit 1;;
esac
