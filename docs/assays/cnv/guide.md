# CNV Segments to Gene Sets (`cnv_gene_extractor`)

This guide covers CNV segment extraction from segment tables to gene-level programs and GMT gene sets.

`cnv_gene_extractor` accepts tumor/germline CNV segments, maps them to genes via gene-body overlap,
applies focality penalties to reduce broad-event dominance, and emits amplification/deletion programs.

## Changelog (CNV extractor)

- Added `--segments_format {auto,gdc_seg,cbio_seg}` with cBio SEG autodetection.
- Added `cnv_program_manifest.tsv` with emitted/skipped program QC stats and tie diagnostics.
- Added `focal` preset and expanded `qc` preset to emit focal+broad variants.
- Added deterministic tie handling diagnostics at the top-k cutoff.

## Quickstart

Single-sample or multi-sample segments table:

```bash
geneset-extractors convert cnv_gene_extractor \
  --segments_tsv path/to/segments.tsv \
  --gtf path/to/genes.gtf.gz \
  --out_dir results/cnv \
  --organism human \
  --genome_build hg38

geneset-extractors validate results/cnv
```

Common table columns are auto-detected (`Chromosome/chrom`, `Start/start`, `End/end`, `Segment_Mean/segmean`).
Common cBioPortal SEG headers (`ID`, `chrom`, `loc.start`, `loc.end`, `seg.mean`) are also auto-detected.

## Resources

CNV extraction currently has no external bundle dependency.

- Required: `--segments_tsv` plus `--gtf` for gene coordinates.
- Optional: `--purity_tsv` for purity-adjusted amplitudes.
- Optional future resources (not required in current defaults): centromere/gap masks, cancer-type reference panels.

This means first runs are fully local as long as segment and GTF files are present.

## Inputs

Required:

- `--segments_tsv`
- `--gtf`
- `--organism`, `--genome_build`

Optional column overrides:

- `--chrom_column`, `--start_column`, `--end_column`, `--amplitude_column`, `--sample_id_column`
- `--segments_format {auto,gdc_seg,cbio_seg}` (default `auto`)
- `--coord_system {one_based_closed,zero_based_half_open}` (default `one_based_closed`)
- `--chrom_prefix_mode {auto,add_chr,drop_chr,none}` (default `auto`)

If sample column is absent, the converter treats input as a single sample (`sample1`).

### Segment format examples

GDC-like segments (auto mode usually works):

```bash
geneset-extractors convert cnv_gene_extractor \
  --segments_tsv path/to/gdc_segments.tsv.gz \
  --segments_format gdc_seg \
  --gtf path/to/genes.gtf.gz \
  --out_dir results/cnv_gdc \
  --organism human \
  --genome_build hg38
```

cBioPortal SEG (no manual column overrides needed):

```bash
geneset-extractors convert cnv_gene_extractor \
  --segments_tsv path/to/data_CNA.seg \
  --segments_format cbio_seg \
  --gtf path/to/genes.gtf.gz \
  --out_dir results/cnv_cbio \
  --organism human \
  --genome_build hg38
```

## Scoring model

Per segment `i`:

- signed amplitude `a_i` from segment mean/log2-ratio
- focal penalty `1 / (1 + (len_i / scale)^alpha)`
- gene-count penalty (`none`, `inv_sqrt`, `inv`)

Per gene `g`:

- overlap fraction `f_ig` from segment/gene-body overlap
- default score aggregates `a_i * focal_penalty_i * gene_count_penalty_i * f_ig` via `weighted_mean`

Directional programs per sample:

- `amp`: `max(score, 0)`
- `del`: `max(-score, 0)`

Default selection/normalization:

- `--select top_k --top_k 200`
- `--normalize within_set_l1`

GMT defaults:

- `--gmt_topk_list 200`
- `--gmt_min_genes 100 --gmt_max_genes 500`
- sets below minimum are skipped unless `--emit_small_gene_sets true`

## Purity correction (optional)

Purity correction is opt-in:

```bash
geneset-extractors convert cnv_gene_extractor \
  --segments_tsv path/to/segments.tsv \
  --gtf path/to/genes.gtf.gz \
  --out_dir results/cnv \
  --organism human \
  --genome_build hg38 \
  --purity_tsv path/to/purity.tsv \
  --use_purity_correction true
```

Boolean flag forms supported:

```bash
# implicit true
--use_purity_correction

# explicit true
--use_purity_correction true

# explicit false
--use_purity_correction false
```

Related flags:

- `--purity_sample_id_column` (default `sample_id`)
- `--purity_value_column` (default `purity`)
- `--purity_floor` (default `0.1`)
- `--max_abs_amplitude` (default `3.0`)

If purity is missing for some samples, the converter warns and skips purity correction for those samples.

## Presets

- `--program_preset default` (default): focal penalties on, emits `amp` + `del`
- `--program_preset focal`: stronger focal penalties and optional max segment-length cap
- `--program_preset broad`: weaker focal penalties and no gene-count penalty (for arm-level/aneuploidy-heavy data)
- `--program_preset qc`: emits both focal and broad variants (`amp_focal`, `del_focal`, `amp_broad`, `del_broad`)
- `--program_preset all`: includes `abs` and supports cohort output when enabled

Explicit program override:

- `--program_methods amp,del`
- `--program_methods amp,del,abs`

## Cohort recurrence sets (optional)

```bash
geneset-extractors convert cnv_gene_extractor \
  --segments_tsv path/to/segments.tsv \
  --gtf path/to/genes.gtf.gz \
  --out_dir results/cnv \
  --organism human \
  --genome_build hg38 \
  --emit_cohort_sets true
```

Controls:

- `--cohort_score_threshold` (default `0.15`)
- `--cohort_min_fraction` (default `0.05`)
- `--cohort_min_samples` (default `5`)

Outputs at `out_dir/cohort/genesets.gmt` when criteria are met.

## Outputs

Grouped output root:

- `manifest.tsv` with one row per sample/program directory
- `manifest.tsv` also records `scored_ok`, `gmt_emitted`, and `emit_reason`
- `cnv_program_manifest.tsv` with emitted + skipped program QC rows and tie diagnostics
  (`status` in `skipped_input|scored_only|suppressed_gmt|emitted_gmt`)
- `sample=<SAMPLE>/program=<PROGRAM>/geneset.tsv`
- `sample=<SAMPLE>/program=<PROGRAM>/geneset.meta.json`
- optional `geneset.full.tsv`, `genesets.gmt`, and run summaries per directory
- optional cohort GMT at `cohort/genesets.gmt`
- `skipped_programs.json` at output root with structured reasons when requested
  program/sample combinations are skipped (for example `no_positive_signal`)

Sign semantics:

- `program=amp`: copy-gain direction
- `program=del`: copy-loss direction

## Interpreting CNV enrichments

CNV-driven enrichments often contain strong positional/cytoband signals. This is expected.

Recommended two-pass interpretation:

1. QC pass: run enrichment with positional/cytoband collections enabled to confirm the CNV
   signal is coherent with known loci.
2. Mechanism pass: run enrichment with positional collections disabled (or deprioritized)
   to focus on pathway and process-level hypotheses.

Use `cnv_program_manifest.tsv` to pick representative samples/programs before deep
interpretation (for example, compare `median_seg_len_used_bp`, `max_abs_amp_used`, and
`n_genes_output`).

## Common pitfalls

- Genome build mismatch between segment coordinates and `--gtf`
- Chromosome prefix mismatch (`chr1` vs `1`) without proper `--chrom_prefix_mode`
- Wrong coordinate convention (`one_based_closed` vs `zero_based_half_open`)
- Thresholds too strict (`--min_abs_amplitude`, `--gmt_min_genes`)
- Broad-event dominance: use `--program_preset broad` when arm-level CNV is expected
- If broad events dominate but focal drivers are expected, try `--program_preset focal`
- If no rows map to GTF chromosomes, CNV preflight fails fast with remediation hints.
- Cytoband/region enrichments are expected QC in CNV runs; interpret them separately from
  mechanism-level pathway interpretation.

## Documentation map

- repository docs index: `docs/geneset-extractors.md`
- compatibility index: `docs/geneset-extractors.md`
- CNV methods note: `docs/assays/cnv/methods.tex`
- methods index: `docs/methods.tex`
