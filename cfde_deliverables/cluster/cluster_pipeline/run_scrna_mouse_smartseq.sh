#!/usr/bin/env bash
# INPUT-SPECIFIED workflow: Allen/BICCN *** MOUSE *** cortex + hippocampus (CTX-HPF), SMART-seq (full-length).
# *** ORGANISM = MOUSE *** — gene programs in MOUSE symbols. Complements the mouse 10x set (full-length, fewer cells).
# Source: Allen public S3 (anonymous, verified 200). Funding: NIH BRAIN Initiative / BICCN.
# matrix.csv ~6.9GB, metadata ~28MB. Cols: sample_name / subclass_label / external_donor_name_label.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
B="https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq"
RM_CSV=1 bash "$HERE/scrna_run_dataset.sh" "allen_biccn_MOUSE_ctx_hpf_smartseq" \
  "$B/matrix.csv" "$B/metadata.csv" \
  "sample_name" "subclass_label" "external_donor_name_label" "counts" "mouse"
