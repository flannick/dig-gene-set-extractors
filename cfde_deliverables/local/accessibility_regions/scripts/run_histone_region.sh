#!/usr/bin/env bash
# Region-level background-corrected accessibility for all histone marks (sequential; each mark = its own
# per-mark background). Keeps peak BEDs. Reuses the pre-downloaded clean refGene (no race).
set -uo pipefail
PY=/Users/gage/Codex/PIGEAN_EAGGL/.venv_kidney_h5/bin/python
REF=/Users/gage/.claude/jobs/32851e29/tmp/refGene_hg38_clean.txt.gz
SCRIPT="$HOME/Claude/proj-valiation-challenge/accessibility_regions/scripts/derive_accessibility_regions.py"
for mk in H3K4me3 H3K27me3 H3K36me3 H3K27ac H3K4me1 H3K9me3 H3K9ac; do
  echo "=== $mk region-level ==="
  ENC_ASSAY="Histone ChIP-seq" ENC_OUTPUT="replicated peaks" ENC_TARGET="$mk" \
    LIB="ENCODE_${mk}_region" REFGENE="$REF" "$PY" "$SCRIPT" || echo "  $mk FAILED (continuing)"
done
echo "HISTONE REGION-LEVEL DONE"
