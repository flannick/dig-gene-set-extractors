#!/bin/bash
# UGER array task: process the N-th fragment file (N = $SGE_TASK_ID) from filelist.txt.
set -euo pipefail
WORK="${WORK:?set WORK}"; OUTDIR="${OUTDIR:?set OUTDIR}"
idx="${SGE_TASK_ID:-${1:?need index}}"
FRAG=$(sed -n "${idx}p" "$WORK/filelist.txt")
[ -z "$FRAG" ] && { echo "no file at index $idx"; exit 0; }
echo "task $idx -> $FRAG"
bash "$(dirname "$0")/catlas_aggregate_sample.sh" "$FRAG"
