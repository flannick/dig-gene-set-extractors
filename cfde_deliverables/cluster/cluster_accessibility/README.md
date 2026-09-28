# Cluster-portable accessibility gene-set scripts

**Pure Python stdlib — runs on your cluster's Python 3.9 as-is. No toolchain, no conda, no deps.**
Only needs: `python3` + internet to `encodeproject.org` + `hgdownload.soe.ucsc.edu`.
**Lean:** each peak/link file is downloaded → parsed → **deleted**; only the small gene sets are kept — so it works where local storage is tight (run it on the cluster).

These produce the SAME contract packages we built locally (`geneset.tsv` + `genesets.gmt` + `geneset.meta.json` + `geneset.provenance.json`), one per biosample, each citing its ENCODE accession (NHGRI, public, derived).

## Scripts
- **`encode_accessibility_batch.py`** — ENCODE peak→gene (promoter TSS±1kb). Does **ATAC-seq** (default) or **DNase-seq** via env. Auto-downloads UCSC refGene hg38 once.
- **`derive_re2g_genesets.py`** — ENCODE-rE2G cCRE→gene → per-biosample **regulatory-target genes** (distal+proximal).
- **`derive_accessible_genes.py`** — single peak file → accessible genes (used by the 4DN MACS path in `cluster_pipeline/`).

## Run examples (cluster)
```bash
mkdir -p run && cd run

# ENCODE ATAC accessible genes (all human biosamples)
OUTDIR=$PWD/atac TMPDIR=$PWD/tmp python3 encode_accessibility_batch.py

# ENCODE DNase accessible genes
ENC_ASSAY="DNase-seq" ENC_OUTPUT="peaks" ENC_LIB="ENCODE_DNase_accessible" \
  OUTDIR=$PWD/dnase TMPDIR=$PWD/tmp python3 encode_accessibility_batch.py

# ENCODE rE2G regulatory-target genes
OUTDIR=$PWD/re2g TMPDIR=$PWD/tmp python3 derive_re2g_genesets.py
```
All are **resume-safe** (re-run to continue; finished biosamples are skipped) and append progress to `<OUTDIR>/batch_log.txt`. Run under `nohup`/`screen` for long batches.

## Compliance
All inputs ENCODE/NHGRI, public, GRCh38. Outputs are **derived + cited** (don't recreate published gene-set products; we derive from public peaks/links). Report any failures (logged per-biosample) and I'll adjust.
