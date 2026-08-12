#!/usr/bin/env bash
# Package a gene-set output dir into the SDG submissions area with an md5 (Linux cluster).
# Usage: package_submission.sh <output_dir> <submission_name>
# Env: SUBMIT_DIR (default the SDG submissions path).
set -uo pipefail
SRC=${1:?output dir to package}; NAME=${2:?submission name}
SUBMIT_DIR=${SUBMIT_DIR:-/humgen/diabetes2/users/ryank/CFDE/geneset_extractors/submissions/SDG}
[ -d "$SRC" ] || { echo "no such output dir: $SRC"; exit 1; }
mkdir -p "$SUBMIT_DIR"
Z="$SUBMIT_DIR/${NAME}.zip"
( cd "$(dirname "$SRC")" && zip -rq "$Z" "$(basename "$SRC")" -x '*/batch_log.txt' '*/pipeline.log' '*/run.log' )
md5sum "$Z" | tee "$Z.md5"
echo "SUBMITTED -> $Z"
