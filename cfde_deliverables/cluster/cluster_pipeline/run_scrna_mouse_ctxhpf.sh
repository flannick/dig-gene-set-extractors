#!/usr/bin/env bash
# INPUT-SPECIFIED workflow: Allen/BICCN *** MOUSE *** cortex + hippocampus (CTX-HPF), 10x snRNA-seq.
# *** ORGANISM = MOUSE (Mus musculus) *** — gene programs will be in MOUSE gene symbols, NOT human.
#     Outputs are named/cited as mouse so they are never confused with the human sets.
# Source: Allen public S3 (anonymous, verified). Funding: NIH BRAIN Initiative / BICCN (NIH by construction).
# SIZE WARNING: matrix.csv ~73GB; CSV->TSV needs transient space, so RM_CSV=1 drops CSVs after conversion
#   (working set ~75GB). Ensure the WORK filesystem has room. Cols: sample_name / subclass_label / external_donor_name_label.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
B="https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x"
RM_CSV=1 bash "$HERE/scrna_run_dataset.sh" "allen_biccn_MOUSE_ctx_hpf_10x" \
  "$B/matrix.csv" "$B/metadata.csv" \
  "sample_name" "subclass_label" "external_donor_name_label" "counts" "mouse"
