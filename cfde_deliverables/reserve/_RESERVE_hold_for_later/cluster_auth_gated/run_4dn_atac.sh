#!/usr/bin/env bash
# 4DN-original ATAC -> accessible genes (CLUSTER / credentialed).
# WHY here not local: 4DN downloads are access-key gated (HTTP 403 anon), AND 4DN ATAC
# processed files are "read positions" (fragments), NOT called peaks -> need MACS peak-calling.
# Prereqs on cluster:
#   - 4DN access key+secret (generate in your 4DN profile): export FDN_KEY=... FDN_SECRET=...
#   - MACS (pip install macs3) in the gsx310 conda env (or any py>=3.7 env)
#   - UCSC refGene hg38 TSS + the peak->gene extractor (accessible_genes/scripts/derive_accessible_genes.py)
# Usage: run_4dn_atac.sh <4DN_FILE_ACC> <biosample> <outdir>
set -uo pipefail
ACC=${1:?4DN file accession (e.g. 4DNFIWOBZ1TH)}; BS=${2:?biosample}; OUT=${3:?outdir}
REFGENE=${REFGENE:?path to refGene_hg38.txt.gz}
EXTRACTOR=${EXTRACTOR:?path to derive_accessible_genes.py}
: "${FDN_KEY:?set FDN_KEY}"; : "${FDN_SECRET:?set FDN_SECRET}"
mkdir -p "$OUT"; READS="$OUT/$ACC.reads.bed.gz"; PEAKS="$OUT/${ACC}_peaks"

echo "[1] download 4DN read positions (access-key auth)"
curl -s -L --user "$FDN_KEY:$FDN_SECRET" \
  "https://data.4dnucleome.org/files-processed/$ACC/@@download/$ACC.bed.gz" -o "$READS" \
  || { echo "download failed (check key/secret + file released)"; exit 1; }

echo "[2] MACS peak-calling (ATAC; Tn5 read positions as BED)"
# macs3 callpeak: -f BED treats each line as a read; ATAC-style params
gunzip -c "$READS" > "$OUT/$ACC.reads.bed"
macs3 callpeak -t "$OUT/$ACC.reads.bed" -f BED -g hs -n "$ACC" --outdir "$PEAKS" \
  --nomodel --shift -75 --extsize 150 -q 0.01 || { echo "macs3 failed"; exit 1; }
rm -f "$OUT/$ACC.reads.bed" "$READS"   # lean: drop fragments after peak-calling

echo "[3] peaks -> accessible genes (reuse our extractor; GRCh38)"
gzip -f "$PEAKS/${ACC}_peaks.narrowPeak"
python "$EXTRACTOR" "$PEAKS/${ACC}_peaks.narrowPeak.gz" "$REFGENE" "$BS" "$ACC" "4DN_$ACC" "$OUT"
echo "DONE: 4DN accessible-gene set for $BS (cite: 4DN $ACC, NIH Common Fund, public; MACS3-called)."

echo "[4] package -> SDG submissions (zip + md5)"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash "$HERE/package_submission.sh" "$OUT" "4dn_atac_accessible_${ACC}"
