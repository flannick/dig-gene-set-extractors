#!/usr/bin/env bash
# ENCODE histone-mark gene sets (4 marks), each per-biosample. Reuses encode_accessibility_batch.py.
# H3K4me3 already done separately (active promoter). Here: repressed / transcribed-body / enhancer marks.
PY=/Users/gage/Codex/PIGEAN_EAGGL/.venv_kidney_h5/bin/python
SCRIPT="$HOME/Claude/proj-valiation-challenge/accessible_genes/scripts/encode_accessibility_batch.py"
REF=/Users/gage/.claude/jobs/32851e29/tmp/refGene_hg38.txt.gz
TMP=/Users/gage/.claude/jobs/32851e29/tmp
BASE="$HOME/Claude/proj-valiation-challenge/histone_genesets"
run(){ local tgt="$1" out="$2" mode="$3" desc="$4" lib="$5"
  echo "=== $tgt (mode=$mode) ==="
  ENC_ASSAY="Histone ChIP-seq" ENC_TARGET="$tgt" ENC_OUTPUT="$out" ENC_MODE="$mode" ENC_DESC="$desc" ENC_LIB="$lib" \
  OUTDIR="$BASE/$tgt" TMPDIR="$TMP" REFGENE="$REF" "$PY" "$SCRIPT"
}
run H3K27me3 "replicated peaks" promoter "with promoter H3K27me3 (Polycomb-repressed)"            ENCODE_H3K27me3_repressed_promoter
run H3K36me3 "replicated peaks" body     "overlapping the gene body with H3K36me3 (transcribed)"   ENCODE_H3K36me3_transcribed_body
run H3K27ac  "replicated peaks" promoter "with promoter-proximal H3K27ac (active; full enhancer reach via rE2G)" ENCODE_H3K27ac_promoter_proximal
run H3K4me1  "replicated peaks" promoter "with promoter-proximal H3K4me1 (enhancer mark; full reach via rE2G)"   ENCODE_H3K4me1_promoter_proximal
echo "ALL HISTONE MARKS DONE"
