# scRNA → cNMF gene-program pipeline (cluster) — SELF-CONTAINED

Produces cNMF gene-program gene sets from public NIH scRNA-seq. **No external toolchain.** The pipeline
is one Python script (`scrna_cnmf_programs.py`) that uses the installed `cnmf` + `anndata` packages and
writes our contract gene sets (`geneset.tsv`/`.gmt`/`meta`/`provenance`) directly.

> Correction: an earlier version of these scripts assumed a `geneset-extractors` CLI with
> `scrna_cnmf_prepare` / `cnmf_select_k` / `rna_sc_programs` subcommands. **That tool does not exist** —
> it was an incorrect inference. Everything now runs on real, installed packages only.

## Locations (this deployment)
- **SH files live in:** `/humgen/diabetes2/users/gage/CFDE` (keep `package_submission.sh` + `scrna_cnmf_programs.py` alongside the run scripts — they're found via the script's own dir).
- **Outputs zip to:** `/humgen/diabetes2/users/ryank/CFDE/geneset_extractors/submissions/SDG` (override with `SUBMIT_DIR=`). Linux `md5sum`.

## Setup (once)
- `cnmf`/`scanpy`/`anndata` need **Python ≥ 3.10** (cluster base is 3.9) → `setup_env.sh` builds a user-space
  env `gsx310` with **micromamba** (low-memory solver) + conda-forge prebuilt stack (no source builds).
- Run `bash setup_env.sh` (needs internet; if a compute node lacks it or the solve is killed for memory,
  run on the **login node** `dig-ae-dev-03`). **Send the `=== 6. VERIFY ===` block** — it confirms the
  scientific imports AND that the `cNMF` methods we use exist on your install.

## Run (input-specified per dataset — verified anon-public + NIH BICCN)
```
bash run_scrna_allen_human_smartseq.sh   # human cortex SMART-seq (~5.4GB)  — smallest, pilot here
bash run_scrna_allen_m1.sh               # human M1 10x (~7.7GB)
bash run_scrna_mouse_ctxhpf.sh           # MOUSE ctx-hpf 10x (~73GB; RM_CSV=1; ORGANISM=mouse)
```
Each downloads matrix+metadata, CSV→TSV, runs cNMF, writes gene sets, and zips them to SDG. To run any
other matrix directly: `bash run_scrna_programs.sh <matrix.tsv> <meta.tsv> <out> <name>` with the column
env vars (`CELL_ID_COLUMN`/`CELL_TYPE_COLUMN`/`DONOR_COLUMN`, plus `ORGANISM`, `K`, `TOPN`, `MAX_CELLS`).

## What the pipeline does (scrna_cnmf_programs.py)
1. read metadata, **balanced subsample** of cells across cell types (`MAX_CELLS`, default 20k);
2. stream the cells×genes matrix, keep selected cells → AnnData (raw counts);
3. **cNMF** `prepare`→`factorize`→`combine`→`consensus` at `K` (default 10), `NHVG` highly-variable genes;
4. take **top-N genes per program** (default 100) from the gene-spectra scores → contract gene sets;
5. package to SDG.
Tunables via env: `K`, `NHVG`, `NITER`, `TOPN`, `MAX_CELLS`, `ORGANISM`, `CITATION`.

## Compliance
All wired datasets are **anon-public + NIH-funded** (Allen/BICCN, NIH BRAIN Initiative; verified HTTP 200).
Provenance on every set cites the dataset + organism + cNMF method. HuBMAP/4DN are auth-gated → reserved,
not used here. (HTAN if ever added: open/Level-3+ only — raw is dbGaP-controlled. GEO: verify NIH funding per study.)

## Honest caveats
- I can't execute on your cluster → the `cnmf` API calls are from the package's documented usage; the
  VERIFY block checks they exist before any long run. Report the pilot output and we iterate.
- `cnmf` factorization is compute-heavy — pilot the SMART-seq set first to gauge runtime before the 73GB mouse set.
