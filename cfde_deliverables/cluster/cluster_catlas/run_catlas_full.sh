#!/bin/bash
# CATLAS full-222 driver.  Fragment-aggregation -> per-cell-type accessibility gene sets.
# Data: GSE184462 (Zhang et al. 2021 Cell, CATlas) — NIH-funded, anonymously public.
#
# Usage:
#   bash run_catlas_full.sh setup     # download RAW.tar(35.6G)+metadata, untar, build promoters, make filelist
#   bash run_catlas_full.sh run       # aggregate all samples (UGER array if available, else local parallel)
#   bash run_catlas_full.sh combine   # combine partials -> gene sets -> package to SDG
#   bash run_catlas_full.sh all       # setup -> run -> combine
#
# Env knobs: WORK (default: this dir), SDG (submission dir), NPROC (local parallelism),
#            MINSIG/MINCOUNT/LOW/HIGH/MIN_CELLS (calling thresholds, see catlas_combine_call.py).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
export WORK="${WORK:-$HERE/catlas_work}"
export OUTDIR="${OUTDIR:-$WORK/partials}"
export META="$WORK/GSE184462_metadata.tsv.gz"
SDG="${SDG:-/humgen/diabetes2/users/ryank/CFDE/geneset_extractors/submissions/SDG}"
NPROC="${NPROC:-8}"
mkdir -p "$WORK" "$OUTDIR"
RAW="$WORK/GSE184462_RAW.tar"
FRAGDIR="$WORK/frags"

dl(){ # url dest
  local url="$1" dest="$2"
  [ -s "$dest" ] && { echo "have $(basename "$dest")"; return; }
  echo "downloading $(basename "$dest") ..."
  wget -q -O "$dest.part" "$url" && mv "$dest.part" "$dest"
}

setup(){
  dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184462/suppl/GSE184462_metadata.tsv.gz" "$META"
  dl "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE184462&format=file" "$RAW"
  echo "RAW.tar size: $(ls -lh "$RAW" | awk '{print $5}')"
  mkdir -p "$FRAGDIR"
  if [ -z "$(ls "$FRAGDIR"/*_fragments.bed.gz 2>/dev/null)" ]; then
    echo "extracting fragment bed.gz from tar ..."
    tar -xf "$RAW" -C "$FRAGDIR" --wildcards '*_fragments.bed.gz'
  fi
  ls "$FRAGDIR"/*_fragments.bed.gz > "$WORK/filelist.txt"
  echo "fragment files: $(wc -l < "$WORK/filelist.txt")"
  echo "building promoter bins ..."
  WORK="$WORK" python3 "$HERE/catlas_make_promoters.py"
}

run(){
  local n; n=$(wc -l < "$WORK/filelist.txt")
  echo "aggregating $n samples ..."
  if command -v qsub >/dev/null 2>&1; then
    echo "submitting UGER array 1-$n (blocking with -sync y) ..."
    qsub -sync y -t 1-"$n" -cwd -j y -o "$WORK/logs/" -l h_vmem=8g \
         -N catlas_agg -v WORK="$WORK",OUTDIR="$OUTDIR" \
         "$HERE/catlas_array_task.sh"
  else
    echo "no qsub; running locally with $NPROC procs ..."
    mkdir -p "$WORK/logs"
    cat "$WORK/filelist.txt" | xargs -P "$NPROC" -I{} bash "$HERE/catlas_aggregate_sample.sh" {}
  fi
  echo "partials written: $(ls "$OUTDIR"/*.partial.tsv 2>/dev/null | wc -l)"
}

combine(){
  echo "combining -> gene sets ..."
  WORK="$WORK" PARTDIR="$OUTDIR" META="$META" \
    OUTDIR_SETS="$WORK/catlas_genesets_full" python3 "$HERE/catlas_combine_call.py"
  # package
  local sets="$WORK/catlas_genesets_full"
  local n; n=$(ls -d "$sets"/*/ 2>/dev/null | wc -l)
  echo "gene sets: $n -> packaging"
  ( cd "$sets/.." && zip -qr "$WORK/batch_catlas_full222.zip" "$(basename "$sets")" )
  if [ -d "$SDG" ]; then cp "$WORK/batch_catlas_full222.zip" "$SDG/"; echo "copied to $SDG"; fi
  md5sum "$WORK/batch_catlas_full222.zip"
}

case "${1:-all}" in
  setup) setup;;
  run) run;;
  combine) combine;;
  all) setup; run; combine;;
  *) echo "usage: $0 {setup|run|combine|all}"; exit 1;;
esac
