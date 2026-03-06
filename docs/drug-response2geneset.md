# Drug Response to Gene Sets (`drug_response_screen`)

This is the practical runbook for drug-response extraction.

Primary CLI:

- `geneset-extractors`
- `omics2geneset` is a compatibility alias in this repo

Methods note: `docs/assays/drug_response/methods.tex`.
Bundle guide: `docs/assays/drug_response/reference_bundle.md`.

## Minimal long-table run

```bash
geneset-extractors convert drug_response_screen \
  --response_tsv path/to/response.tsv \
  --drug_targets_tsv path/to/drug_targets.tsv \
  --groups_tsv path/to/groups.tsv \
  --organism human \
  --genome_build hg38 \
  --out_dir out/drug_response
```

Required flags in practice:

- `--organism`
- `--genome_build`
- `--out_dir`
- one input mode:
  - generic long tables (`--response_tsv` + `--drug_targets_tsv`), or
  - PRISM mode (`--prism_matrix_csv` + `--prism_treatment_info_csv` + `--prism_cell_line_info_csv`)

## Input contracts

### `--response_tsv` (long format)

Required columns:

- `sample_id` (string)
- `drug_id` (string)
- `response` (float)

### `--drug_targets_tsv`

Required columns:

- `drug_id`
- `gene_symbol`

Optional:

- `weight`
- `source`

### `--groups_tsv`

Required columns:

- `sample_id`
- `group`

## PRISM preparation workflow

Use one command to fetch and standardize PRISM inputs into long-form tables:

```bash
geneset-extractors workflows prism_prepare \
  --out_dir out/prism_prepare
```

Default file-id template is Figshare (`https://ndownloader.figshare.com/files/{file_id}`).
If a template returns JSON/HTML metadata instead of tabular data, workflow
summary and warnings make this explicit and suggest the Figshare template.

This writes:

- `response_long.tsv` with columns `sample_id`, `drug_id`, `response`
- `groups.tsv` with columns `sample_id`, `group` (`primary_tissue`)
- `drug_targets.tsv` with columns `drug_id`, `gene_symbol`, `weight`
- `prepare_summary.json` with fetch and parse details

Quick deterministic subset for rapid checks:

```bash
geneset-extractors workflows prism_prepare \
  --out_dir out/prism_prepare_quick \
  --subset_seed 11 \
  --balance_by primary_tissue \
  --min_per_balance_bin 5 \
  --max_cell_lines_per_group 50 \
  --max_cell_lines_total 500 \
  --max_compounds_total 300
```

Then run extraction:

```bash
geneset-extractors convert drug_response_screen \
  --response_tsv out/prism_prepare/response_long.tsv \
  --drug_targets_tsv out/prism_prepare/drug_targets.tsv \
  --groups_tsv out/prism_prepare/groups.tsv \
  --organism human \
  --genome_build hg38 \
  --out_dir out/prism_programs
```

With an optional local bundle for more connectable defaults:

```bash
geneset-extractors convert drug_response_screen \
  --response_tsv out/prism_prepare/response_long.tsv \
  --groups_tsv out/prism_prepare/groups.tsv \
  --organism human \
  --genome_build hg38 \
  --resources_dir /path/to/drug_response_bundle \
  --out_dir out/prism_programs
```

## Defaults that reduce spurious outputs

- `--min_group_size 5`:
  groups smaller than this are skipped with explicit warnings and metadata.
- `--max_targets_per_drug 50 --target_promiscuity_policy warn`:
  flags very large target lists; switch to `drop` or `cap` if needed.
- `--target_ubiquity_penalty idf`:
  downweights genes targeted by many drugs.
- response-side ubiquity:
  - preferred flag: `--response_ubiquity_penalty`
  - legacy-compatible alias: `--ubiquity_penalty`

## Output and QC summary

Root outputs include:

- `manifest.tsv` (emitted program rows only)
- `group_qc.tsv` (all groups with `sets_emitted` and `reason_if_skipped`)
- `run_summary.json`
- `run_summary.txt`
- optional combined `genesets.gmt`

Per-program outputs include:

- `program=<ID>/geneset.tsv`
- `program=<ID>/geneset.meta.json`
- optional `program=<ID>/genesets.gmt`

The run summary reports:

- sample/drug/response counts
- groups total/emitted/skipped and per-group counts
- target coverage and promiscuity stats
- emitted set counts and size distribution

## GMT compatibility mode

Default format is DIG 2-column:

- `set_id<TAB>gene1 gene2 ...`

For classic GMT parsers use:

```bash
--gmt_format classic
```

Classic format:

- `set_id<TAB>na<TAB>gene1<TAB>gene2...`

## Common failure modes

- Missing columns:
  check `geneset-extractors convert drug_response_screen --help` and align column flags.
- Tiny groups:
  if a group is skipped, either increase group size or lower `--min_group_size`.
- PRISM download produced JSON/HTML:
  check `prepare_summary.json` fetch `sniff` fields; use
  `--depmap_file_id_url_template https://ndownloader.figshare.com/files/{file_id}`.
- Target table issues:
  missing targets or overly large target lists reduce specificity.
  Inspect `run_summary.json` target stats and adjust promiscuity settings.
- Broad GPCR/neuroactive enrichments:
  common in pan-pharmacology settings; tighten target filtering and keep IDF penalty enabled.
