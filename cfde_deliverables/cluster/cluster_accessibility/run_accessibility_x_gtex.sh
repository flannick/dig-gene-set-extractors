#!/usr/bin/env bash
# CLUSTER: accessibility (ATAC/DNase) consensus x GTEx tissue-enriched, against ALREADY-UPLOADED sets
# (no re-download). Then package into SDG submissions. Pure-stdlib python (3.9+ ok).
# Edit the two *_DIR paths to point at the uploaded set directories on the cluster.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY=${PY:-python3}

# ---- EDIT THESE to your cluster locations ----
ATAC_DIR=${ATAC_DIR:-/humgen/diabetes2/users/gage/CFDE/sets/ENCODE_ATAC_accessible}     # per-biosample ATAC set dirs
DNASE_DIR=${DNASE_DIR:-/humgen/diabetes2/users/gage/CFDE/sets/ENCODE_DNase_accessible}  # per-biosample DNase set dirs
ENRICHED_DIR=${ENRICHED_DIR:-/humgen/diabetes2/users/gage/CFDE/sets/GTEx_tissue_enriched} # GTEx tissue-enriched set dirs
OUTBASE=${OUTBASE:-$HERE/output}
# ----------------------------------------------

run(){ # <access_dir> <libname>
  local ad="$1" lib="$2" out="$OUTBASE/$lib"
  [ -d "$ad" ] || { echo "SKIP $lib: missing $ad"; return; }
  ACCESS_DIR="$ad" ENRICHED_DIR="$ENRICHED_DIR" OUTDIR="$out" LIBNAME="$lib" \
    "$PY" "$HERE/derive_accessibility_x_gtex.py"
  bash "$HERE/package_submission.sh" "$out" "${lib}_x_GTEx"
}
run "$ATAC_DIR"  "ENCODE_ATAC_accessible"
run "$DNASE_DIR" "ENCODE_DNase_accessible"
echo "accessibility x GTEx DONE"
