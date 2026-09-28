#!/usr/bin/env bash
# Additional ENCODE histone marks: H3K9me3 (constitutive heterochromatin -> silenced) and H3K9ac (active).
PY=/Users/gage/Codex/PIGEAN_EAGGL/.venv_kidney_h5/bin/python
SCRIPT="$HOME/Claude/proj-valiation-challenge/accessible_genes/scripts/encode_accessibility_batch.py"
REF=/Users/gage/.claude/jobs/32851e29/tmp/refGene_hg38.txt.gz
TMP=/Users/gage/.claude/jobs/32851e29/tmp
BASE="$HOME/Claude/proj-valiation-challenge/histone_genesets"
run(){ ENC_ASSAY="Histone ChIP-seq" ENC_TARGET="$1" ENC_OUTPUT="replicated peaks" ENC_MODE=promoter ENC_DESC="$2" ENC_LIB="$3" \
       OUTDIR="$BASE/$1" TMPDIR="$TMP" REFGENE="$REF" "$PY" "$SCRIPT"; }
run H3K9me3 "with promoter H3K9me3 (constitutive heterochromatin / silenced)" ENCODE_H3K9me3_silenced_promoter
run H3K9ac  "with promoter-proximal H3K9ac (active)"                          ENCODE_H3K9ac_active_promoter
echo "H3K9 MARKS DONE"
