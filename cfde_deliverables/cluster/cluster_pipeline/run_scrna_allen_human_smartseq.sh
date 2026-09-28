#!/usr/bin/env bash
# INPUT-SPECIFIED workflow: Allen/BICCN human cortex, SMART-seq (full-length) snRNA-seq.
# Source: Allen public S3 (anonymous, verified). Funding: NIH BRAIN Initiative / BICCN (NIH by construction).
# matrix.csv ~5.4GB (cells x genes, read counts), metadata.csv ~14MB. Cols: sample_name / subclass_label / external_donor_name_label.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
B="https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq"
bash "$HERE/scrna_run_dataset.sh" "allen_biccn_human_ctx_smartseq" \
  "$B/matrix.csv" "$B/metadata.csv" \
  "sample_name" "subclass_label" "external_donor_name_label" "counts"
