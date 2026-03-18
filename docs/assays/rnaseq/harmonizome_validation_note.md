# Harmonizome DE Preset Validation

Date: 2026-03-17

Checkout base: `56881f3468f65c8a71f7fa8d98ea8eb12045f72f`

Scope:
- validate the workflow-side `--de_mode harmonizome` preset after allowing explicit covariates
- compare a real GTEx adipose contrast in:
  - general all-sample mode with covariates
  - conservative balanced Harmonizome mode with covariates
- confirm that the conservative preset is explicit, auditable, and not dominated by the prior mitochondrial / housekeeping failure mode
- run lighter sanity checks on GTEx brain and blood

Important note:
- this note is about DE preparation, not extractor-side post-processing
- the comparison below is `modern + covariates` versus `harmonizome + covariates`
- earlier no-covariate runs were materially more generic; the point of this pass was to validate the new conservative preset that combines balancing with explicit covariates

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
- `78 passed in 15.95s`

Real adipose all-sample run with explicit covariates:

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
  --covariates SEX,SMTSD \
  --backend r_limma_voom \
  --de_mode modern \
  --out_dir tests/tmp/harmonizome_validation/adipose_modern_covariates \
  --organism human \
  --genome_build hg38
```

Real adipose conservative Harmonizome run:

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
  --covariates SEX,SMTSD \
  --backend r_limma_voom \
  --de_mode harmonizome \
  --balance_seed 1 \
  --out_dir tests/tmp/harmonizome_validation/adipose_harmonizome_covariates \
  --organism human \
  --genome_build hg38
```

Light brain and blood sanity checks:

```bash
cd /Users/flannick/codex-workspace/code/dig_gene_set_extractors_module/dig-gene-set-extractors
R_LIBS_USER=/Users/flannick/codex-workspace/code/.Rlib \
PYTHONPATH=src ../../.venv/bin/python -m geneset_extractors.cli workflows rna_de_prepare \
  --modality bulk \
  --counts_tsv tests/tmp/harmonizome_validation/Brain.tsv \
  --matrix_orientation gene_by_sample \
  --feature_id_column gene_id \
  --matrix_gene_symbol_column gene_symbol \
  --sample_metadata_tsv /Users/flannick/codex-workspace/analysis/gene_set_extractors/rna_seq_gene_extractor/prep/gtex_sample_metadata.tsv \
  --sample_id_column SAMPID \
  --subject_metadata_tsv /Users/flannick/codex-workspace/analysis/gene_set_extractors/rna_seq_gene_extractor/prep/gtex_subject_metadata.tsv \
  --subject_join_sample_column SUBJID \
  --subject_join_metadata_column SUBJID \
  --comparisons_tsv tests/tmp/harmonizome_validation/gtex_brain_20_29_vs_60_69.tsv \
  --stratify_by SMTS \
  --covariates SEX,SMTSD \
  --backend r_limma_voom \
  --de_mode harmonizome \
  --balance_seed 1 \
  --out_dir tests/tmp/harmonizome_validation/brain_harmonizome_covariates \
  --organism human \
  --genome_build hg38

R_LIBS_USER=/Users/flannick/codex-workspace/code/.Rlib \
PYTHONPATH=src ../../.venv/bin/python -m geneset_extractors.cli workflows rna_de_prepare \
  --modality bulk \
  --counts_tsv tests/tmp/harmonizome_validation/Blood.tsv \
  --matrix_orientation gene_by_sample \
  --feature_id_column gene_id \
  --matrix_gene_symbol_column gene_symbol \
  --sample_metadata_tsv /Users/flannick/codex-workspace/analysis/gene_set_extractors/rna_seq_gene_extractor/prep/gtex_sample_metadata.tsv \
  --sample_id_column SAMPID \
  --subject_metadata_tsv /Users/flannick/codex-workspace/analysis/gene_set_extractors/rna_seq_gene_extractor/prep/gtex_subject_metadata.tsv \
  --subject_join_sample_column SUBJID \
  --subject_join_metadata_column SUBJID \
  --comparisons_tsv tests/tmp/harmonizome_validation/gtex_blood_20_29_vs_60_69.tsv \
  --stratify_by SMTS \
  --covariates SEX,SMTSD \
  --backend r_limma_voom \
  --de_mode harmonizome \
  --balance_seed 1 \
  --out_dir tests/tmp/harmonizome_validation/blood_harmonizome_covariates \
  --organism human \
  --genome_build hg38
```

## Key Results

### Real contrast: `GTEx_AdiposeTissue_20-29_vs_50-59`

From `comparison_audit.tsv` and `deg_long.tsv`:

| mode | covariates | pre-balance counts | post-balance counts | selected samples | tested genes | significant genes (`padj <= 0.05`) |
|---|---|---:|---:|---:|---:|---:|
| `modern` | `SEX,SMTSD` | `397 vs 97` | `397 vs 97` | `494` | `23,426` | `12,346` |
| `harmonizome` | `SEX,SMTSD` | `397 vs 97` | `97 vs 97` | `194` | `20,999` | `9,376` |

Resolved preset audit details:
- `covariates_used = SEX,SMTSD` in both runs
- `harmonizome_covariate_mode = explicit` in the balanced run
- `balance_seed = 1` in the balanced run
- `comparison_selected_samples.tsv` recorded all `194` chosen sample IDs for the balanced run

Top positive genes by adjusted significance:

| mode | top positive genes |
|---|---|
| `modern` + covariates | `EYA4`, `PYHIN1`, `EDA2R`, `PTCHD4`, `ZMAT3`, `SMOC2`, `LMO3`, `CXorf57`, `TNFRSF14`, `LINC01759` |
| `harmonizome` + covariates | `EDA2R`, `PYHIN1`, `EYA4`, `ZMAT3`, `PTCHD4`, `SMOC2`, `GAS6`, `ABCB5`, `LINC01759`, `GALNT8` |

Specific observation:
- neither top-10 list is dominated by `MT-*`, `EEF1A1`, `ACTB`, `ACTG1`, or `TMSB4X`
- that is a clear improvement over the earlier imbalanced/no-covariate failure mode

Interpretation:
- adding explicit covariates already improves the general all-sample fit materially
- the balanced Harmonizome preset is still the more conservative setting:
  - fewer retained samples
  - fewer tested genes
  - fewer significant genes
- this is the intended behavior: the preset is designed for conservative signature generation in broad tissues, not for maximum-power DE

### Light sanity checks

Balanced Harmonizome runs with `SEX,SMTSD` for two additional tissues:

| contrast | pre-balance counts | post-balance counts | tested genes | significant genes (`padj <= 0.05`) | top positive genes |
|---|---:|---:|---:|---:|---|
| `GTEx_Brain_20-29_vs_60-69` | `1317 vs 70` | `70 vs 70` | `21,661` | `9,919` | `RP11-317N12.1`, `ATP6AP1L`, `HHLA3`, `GFAP`, `CSTB`, `SGF29`, `ITGB4`, `MACROD1`, `GSTP1`, `MT1G` |
| `GTEx_Blood_20-29_vs_60-69` | `295 vs 85` | `85 vs 85` | `18,045` | `8,709` | `FCRL6`, `B3GAT1`, `S1PR5`, `GZMB`, `PRSS23`, `LGR6`, `GZMH`, `ADGRG1`, `PPP2R2B`, `PYHIN1` |

Specific observation:
- these top lists are not dominated by the earlier generic housekeeping set (`EEF1A1`, `ACTB`, `ACTG1`, `TMSB4X`)
- brain includes one metal-response gene (`MT1G`) but not a broad mitochondrial takeover
- blood is NK / immune flavored rather than obviously housekeeping-dominated

## Conclusion

The conservative preset is now behaving as intended:

- it is explicit rather than hidden
- it balances each contrast deterministically
- it records the selected samples, counts, covariates, and seed
- it supports the recommended broad-tissue covariates (`SEX,SMTSD`)
- it remains distinct from the general `modern` mode

Most importantly, the adipose validation no longer looks like the earlier generic mitochondrial / housekeeping failure. The practical effect of the new preset is not that it magically creates tissue specificity from nothing; it is that it provides a narrower, more conservative DE fit for signature construction when broad-tissue heterogeneity would otherwise dominate.
