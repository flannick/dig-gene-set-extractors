# Harmonizome DE Mode Validation

Date: 2026-03-17

Checkout base: `df7a285cf1dd4e044f93494ae4d8bdfdb4d7abd9`

Scope:
- validate the new workflow-side `--de_mode harmonizome` preset
- confirm deterministic per-comparison balancing
- compare a real GTEx adipose contrast before vs after balancing
- confirm the multi-contrast GTEx-style path still executes

Important note:
- the raw notebook DE table for `GTEx_AdiposeTissue_20-29_vs_50-59` was not available locally in this checkout
- the public downloaded Harmonizome GMT was used as the public reference artifact for overlap checks

## Commands Run

Focused RNA regressions:

```bash
cd /Users/flannick/codex-workspace/code/dig_gene_set_extractors_module/dig-gene-set-extractors
PYTHONPATH=src ../../.venv/bin/python -m pytest \
  tests/test_rna_de_prepare_workflow.py \
  tests/test_rna_deg_converter.py \
  tests/test_rna_deg_multi_converter.py \
  tests/test_cli_smoke.py \
  -q
```

Result:
- `76 passed in 16.12s`

Real single-contrast modern run:

```bash
cd /Users/flannick/codex-workspace/code/dig_gene_set_extractors_module/dig-gene-set-extractors
R_LIBS_USER=/Users/flannick/codex-workspace/code/.Rlib \
PYTHONPATH=src ../../.venv/bin/python -m geneset_extractors.cli workflows rna_de_prepare \
  --modality bulk \
  --counts_tsv tests/tmp/harmonizome_validation/Adipose_Tissue.tsv \
  --matrix_orientation gene_by_sample \
  --feature_id_column gene_id \
  --matrix_gene_symbol_column gene_symbol \
  --sample_metadata_tsv /Users/flannick/codex-workspace/analysis/gene_set_extractors/rna_seq_gene_extractor/prep/gtex_sample_metadata.tsv \
  --sample_id_column SAMPID \
  --subject_metadata_tsv /Users/flannick/codex-workspace/analysis/gene_set_extractors/rna_seq_gene_extractor/prep/gtex_subject_metadata.tsv \
  --subject_join_sample_column SUBJID \
  --subject_join_metadata_column SUBJID \
  --comparisons_tsv tests/tmp/harmonizome_validation/gtex_adipose_20_29_vs_50_59.tsv \
  --stratify_by SMTS \
  --backend r_limma_voom \
  --de_mode modern \
  --out_dir tests/tmp/harmonizome_validation/adipose_modern \
  --organism human \
  --genome_build hg38
```

Real single-contrast Harmonizome run:

```bash
cd /Users/flannick/codex-workspace/code/dig_gene_set_extractors_module/dig-gene-set-extractors
R_LIBS_USER=/Users/flannick/codex-workspace/code/.Rlib \
PYTHONPATH=src ../../.venv/bin/python -m geneset_extractors.cli workflows rna_de_prepare \
  --modality bulk \
  --counts_tsv tests/tmp/harmonizome_validation/Adipose_Tissue.tsv \
  --matrix_orientation gene_by_sample \
  --feature_id_column gene_id \
  --matrix_gene_symbol_column gene_symbol \
  --sample_metadata_tsv /Users/flannick/codex-workspace/analysis/gene_set_extractors/rna_seq_gene_extractor/prep/gtex_sample_metadata.tsv \
  --sample_id_column SAMPID \
  --subject_metadata_tsv /Users/flannick/codex-workspace/analysis/gene_set_extractors/rna_seq_gene_extractor/prep/gtex_subject_metadata.tsv \
  --subject_join_sample_column SUBJID \
  --subject_join_metadata_column SUBJID \
  --comparisons_tsv tests/tmp/harmonizome_validation/gtex_adipose_20_29_vs_50_59.tsv \
  --stratify_by SMTS \
  --backend r_limma_voom \
  --de_mode harmonizome \
  --out_dir tests/tmp/harmonizome_validation/adipose_harmonizome \
  --organism human \
  --genome_build hg38
```

Post-process both DE tables with the current extractor defaults:

```bash
cd /Users/flannick/codex-workspace/code/dig_gene_set_extractors_module/dig-gene-set-extractors
PYTHONPATH=src ../../.venv/bin/python -m geneset_extractors.cli convert rna_deg_multi \
  --deg_tsv tests/tmp/harmonizome_validation/adipose_modern/deg_long.tsv \
  --comparison_column comparison_id \
  --out_dir tests/tmp/harmonizome_validation/adipose_modern_extract \
  --organism human \
  --genome_build hg38

PYTHONPATH=src ../../.venv/bin/python -m geneset_extractors.cli convert rna_deg_multi \
  --deg_tsv tests/tmp/harmonizome_validation/adipose_harmonizome/deg_long.tsv \
  --comparison_column comparison_id \
  --out_dir tests/tmp/harmonizome_validation/adipose_harmonizome_extract \
  --organism human \
  --genome_build hg38
```

GTEx-style multi-contrast smoke with the full age-comparisons table:

```bash
cd /Users/flannick/codex-workspace/code/dig_gene_set_extractors_module/dig-gene-set-extractors
PYTHONPATH=src ../../.venv/bin/python -m geneset_extractors.cli workflows rna_de_prepare \
  --modality bulk \
  --counts_tsv tests/tmp/harmonizome_validation/Adipose_Tissue.tsv \
  --matrix_orientation gene_by_sample \
  --feature_id_column gene_id \
  --matrix_gene_symbol_column gene_symbol \
  --sample_metadata_tsv /Users/flannick/codex-workspace/analysis/gene_set_extractors/rna_seq_gene_extractor/prep/gtex_sample_metadata.tsv \
  --sample_id_column SAMPID \
  --subject_metadata_tsv /Users/flannick/codex-workspace/analysis/gene_set_extractors/rna_seq_gene_extractor/prep/gtex_subject_metadata.tsv \
  --subject_join_sample_column SUBJID \
  --subject_join_metadata_column SUBJID \
  --comparisons_tsv /Users/flannick/codex-workspace/analysis/gene_set_extractors/rna_seq_gene_extractor/prep/gtex_age_comparisons.tsv \
  --stratify_by SMTS \
  --backend lightweight \
  --de_mode modern \
  --out_dir tests/tmp/harmonizome_validation/adipose_full_table_smoke \
  --organism human \
  --genome_build hg38
```

## Key Results

### Real contrast: `GTEx_AdiposeTissue_20-29_vs_50-59`

From `comparison_audit.tsv` and `deg_long.tsv`:

| mode | backend | pre-balance counts | post-balance counts | selected samples | tested genes | significant genes (`padj <= 0.05`) |
|---|---|---:|---:|---:|---:|---:|
| `modern` | `r_limma_voom` | `397 vs 97` | `397 vs 97` | `494` | `23,426` | `11,656` |
| `harmonizome` | `r_limma_voom` | `397 vs 97` | `97 vs 97` | `194` | `20,999` | `8,584` |

Additional balanced-run details:
- `balance_seed = 1`
- `comparison_selected_samples.tsv` recorded all `194` chosen sample IDs
- each group contributed exactly `97` samples after balancing

Relative change from enabling the preset:
- significant genes reduced by `3,072` (`11,656 -> 8,584`)
- tested genes reduced by `2,427` (`23,426 -> 20,999`)

Interpretation:
- this is a clear upstream calibration shift in the notebook direction
- it does not fully close the gap to the reported notebook figure (`6,915` significant genes), so balancing is important but not sufficient on its own

### Public Harmonizome reference overlap

Public reference artifact:
- `/Users/flannick/codex-workspace/analysis/resources/pigean/data/GTEx_XMT_2022-06-06_GTEx_Aging_Signatures_2021.gmt.gz`

Using the top `250` positive and negative genes from `deg_long.tsv`:

| mode | top-250 Up overlap vs public GMT | top-250 Down overlap vs public GMT |
|---|---:|---:|
| `modern` | `27 / 250` | `33 / 250` |
| `harmonizome` | `28 / 250` | `27 / 250` |

Using the top `50` positive and negative genes from `deg_long.tsv`:

| mode | top-50 Up overlap vs public GMT | top-50 Down overlap vs public GMT |
|---|---:|---:|
| `modern` | `6 / 50` | `8 / 50` |
| `harmonizome` | `12 / 50` | `9 / 50` |

Representative top-50 Harmonizome-mode overlaps with the public `Up` set:
- `PYHIN1`
- `EDA2R`
- `EYA4`
- `PTCHD4`
- `GALNT8`
- `LMO3`
- `HOXD9`
- `TRAT1`
- `TSC22D1-AS1`
- `GDF15`
- `EOMES`
- `APOBEC3H`

Representative top-50 Harmonizome-mode overlaps with the public `Down` set:
- `TUBB6`
- `CSPG5`
- `TUBB4B`
- `FNDC4`
- `FGFRL1`
- `ANKRD53`
- `GGCT`
- `PLCXD1`
- `DUSP4`

The extractor post-processing step was transparent here:
- converting `deg_long.tsv` with the current `rna_deg_multi` defaults preserved the same top-250 overlap counts as the raw ranked DE table
- this means the workflow-side DE fit is the main driver of the observed change, not DEG post-processing drift

## Multi-contrast GTEx-style Smoke

The full GTEx age-comparisons table was exercised against the adipose split matrix:
- comparison table rows seen: `135`
- emitted adipose contrasts: `5`
- skipped non-adipose contrasts: `130`
- skip reason: `insufficient_units_after_selection`

This does not constitute a full 30-tissue GTEx rerun, but it does confirm that the new workflow changes did not break the broader multi-contrast GTEx-style path.

## Validation-time fixes discovered and addressed

Two real issues only surfaced under the real GTEx matrix:

1. R backend count staging bug
- the backend input matrix still included the `gene_symbol` column
- `edgeR::DGEList` then saw a mixed metadata/count matrix and failed
- fixed by teaching both `r_limma_voom` and `r_dream` backends to drop `gene_symbol` from the numeric count matrix while preserving it for output

2. CLI consistency bug
- `--group_column` was still required even when every row in `--comparisons_tsv` already provided `group_column`
- fixed so explicit comparisons TSV rows can supply `group_column` without a redundant workflow-level flag

## Bottom Line

The new `--de_mode harmonizome` preset is working as intended:
- explicit
- deterministic
- auditable
- isolated to the workflow layer
- not the universal default

On the real adipose aging contrast, enabling the preset:
- matched the requested balanced sample pool exactly (`97 vs 97`)
- materially reduced the number of significant genes
- moved the fit closer to the published notebook behavior
- improved top-50 positive overlap against the public Harmonizome signature artifact

The improvement is real, but it is partial:
- the balanced run is still above the notebook’s reported significant-gene count
- simple public-GMT overlap does not show a uniformly large gain in every direction

Current conclusion:
- balancing was a necessary upstream correction
- it reduces the generic inflation seen in the unbalanced pool
- additional notebook-parity work likely still lives in the DE fitting details beyond balancing alone
