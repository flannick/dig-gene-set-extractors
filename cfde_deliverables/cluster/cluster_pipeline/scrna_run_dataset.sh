#!/usr/bin/env bash
# Helper: fetch a PUBLIC scRNA matrix+metadata, CSV->TSV, run the cNMF gene-program pipeline, package->SDG.
# Args: NAME MATRIX_URL META_URL CELL_ID_COL CELL_TYPE_COL DONOR_COL VALUE_TYPE ORGANISM
# Matrix is converted with tr (numeric, no quoting); metadata via python csv (handles quoted commas).
set -uo pipefail
NAME=${1:?name}; MURL=${2:?matrix_url}; EURL=${3:?meta_url}
CID=${4:?cell_id_col}; CT=${5:?cell_type_col}; DN=${6:?donor_col}; VT=${7:-counts}; ORG=${8:-human}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK=${WORK:-$HERE/work/$NAME}; OUT=${OUT:-$HERE/out/$NAME}; mkdir -p "$WORK"

echo "[$NAME 1/3] download (anonymous, public)"
# ATOMIC: download to .part then mv — so a file only exists when COMPLETE (prevents a concurrent/2nd
# run from reading a half-written matrix and parsing a truncated cell set, which corrupted M1).
[ -s "$WORK/matrix.csv" ] || { curl -fSL "$MURL" -o "$WORK/matrix.csv.part" && mv "$WORK/matrix.csv.part" "$WORK/matrix.csv"; } || { echo "matrix download failed"; exit 1; }
[ -s "$WORK/meta.csv" ]   || { curl -fSL "$EURL" -o "$WORK/meta.csv.part"   && mv "$WORK/meta.csv.part"   "$WORK/meta.csv"; }   || { echo "meta download failed"; exit 1; }

echo "[$NAME 2/3] CSV -> TSV (atomic)"
[ -s "$WORK/matrix.tsv" ] || { tr ',' '\t' < "$WORK/matrix.csv" > "$WORK/matrix.tsv.part" && mv "$WORK/matrix.tsv.part" "$WORK/matrix.tsv"; }
[ -s "$WORK/meta.tsv" ]   || { python3 -c "import csv,sys;w=csv.writer(sys.stdout,delimiter='\t');[w.writerow(r) for r in csv.reader(open(sys.argv[1]))]" "$WORK/meta.csv" > "$WORK/meta.tsv.part" && mv "$WORK/meta.tsv.part" "$WORK/meta.tsv"; }
# RM_CSV=1 deletes the source CSVs once TSVs exist (saves disk for very large matrices, e.g. the 73GB mouse set)
if [ "${RM_CSV:-0}" = "1" ] && [ -s "$WORK/matrix.tsv" ] && [ -s "$WORK/meta.tsv" ]; then
  rm -f "$WORK/matrix.csv" "$WORK/meta.csv"; echo "  (removed source CSVs; RM_CSV=1)"
fi

echo "[$NAME 3/3] cNMF gene-program pipeline (-> packages to SDG)"
CELL_ID_COLUMN="$CID" CELL_TYPE_COLUMN="$CT" DONOR_COLUMN="$DN" \
ORGANISM="$ORG" CITATION="$NAME (Allen/BICCN; NIH BRAIN Initiative; public)" \
  bash "$HERE/run_scrna_programs.sh" "$WORK/matrix.tsv" "$WORK/meta.tsv" "$OUT" "$NAME"
