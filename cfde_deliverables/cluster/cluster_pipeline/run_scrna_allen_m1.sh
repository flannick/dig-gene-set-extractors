#!/usr/bin/env bash
# INPUT-SPECIFIED workflow: Allen/BICCN human primary motor cortex (M1), 10x snRNA-seq.
# Source: Allen public S3 (anonymous, verified). Funding: NIH BRAIN Initiative / BICCN (NIH by construction).
# matrix.csv ~7.7GB (cells x genes, raw counts), metadata.csv ~24MB. Cols: sample_name / subclass_label / external_donor_name_label.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
B="https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x"
bash "$HERE/scrna_run_dataset.sh" "allen_biccn_human_m1_10x" \
  "$B/matrix.csv" "$B/metadata.csv" \
  "sample_name" "subclass_label" "external_donor_name_label" "counts"
