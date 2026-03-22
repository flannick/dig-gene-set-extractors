# Provenance Compliance Execution Plan

## Goal

Bring all active `dig-gene-set-extractors` extractors onto the current shared provenance contract required by the graph-based REVEAL UI without changing outward artifact names.

## Contract to enforce

- Single-output extractors emit `geneset.tsv`, `geneset.meta.json`, and `geneset.provenance.json`.
- Grouped outputs keep `manifest.tsv` and include `path`, `geneset_id`, `label`, `meta_path`, `provenance_path`, and `focus_node_id`.
- Overlay support flows through `--provenance_overlay_json`.
- Runtime provenance is activated with `activate_runtime_context(...)`.
- Primary inputs use `input_file_record(...)`.
- Bundle/reference resources use `input_file_record(..., resource_record=record)` whenever structured resource metadata exists.
- Grouped manifests are enriched via `enrich_manifest_row(...)`.
- Validation uses the repo’s shared validator and existing metadata writer behavior.

## Work sequence

### Phase 1: Runtime-context gaps

Patch extractors that currently write metadata without activating runtime provenance:

- `atac_bulk`
- `atac_bulk_matrix`
- `proteomics_diff`
- `splice_event_diff`
- `calr_ontology_mapper`
- `calr_profile_query`
- `methylation_dmr`
- `cnv_gene_extractor`
- `morphology_profile_query`
- `chipseq_peak`

Exit criteria:
- each target activates runtime provenance at the top-level entry point that writes metadata or drives metadata construction

### Phase 2: Resource-record enrichment

Patch extractors that already have some provenance support but still collapse bundle/reference resources to local-path-only nodes:

- `atac_bulk`
- `atac_bulk_matrix`
- `atac_sc_10x`
- `splice_event_diff`
- `splice_event_matrix`
- `drug_response_screen`
- `calr_ontology_mapper`
- `calr_profile_query`
- `morphology_profile_query`

Exit criteria:
- resource-manager-backed files are recorded with `resource_record=...`
- provenance file/resource nodes can expose public metadata rather than only local paths

### Phase 3: Manifest and overlay regression coverage

Add or refresh smoke/regression tests for:

- `rna_deg`
- `rna_deg_multi`
- `rna_sc_programs`
- `sc_rna_marker`
- `ptm_site_diff`
- `ptm_site_matrix`
- `methylation_cpg_diff`
- `methylation_dmr_regions`

Exit criteria:
- grouped manifest tests assert enriched columns where applicable
- overlay propagation is covered where requested
- at least one per-extractor output directory validates successfully

### Phase 4: Validation sweep

Run the targeted smoke suite and a shared validation sweep across produced outputs.

Exit criteria:
- patched extractors pass their smoke tests
- validator confirms the expected output contract

## Current status

- Completed runtime-context activation patches for `atac_bulk`, `atac_bulk_matrix`, `proteomics_diff`, `chipseq_peak`, `cnv_gene_extractor`, `methylation_dmr`, `calr_ontology_mapper`, `calr_profile_query`, and `morphology_profile_query`.
- Completed resource-record enrichment patches for `atac_bulk`, `atac_bulk_matrix`, `atac_sc_10x`, `splice_event_diff`, `splice_event_matrix`, `drug_response_screen`, the calorimetry converters/workflows, and `morphology_profile_query`.
- Added or refreshed provenance-focused smoke/regression coverage across ATAC, splice, RNA, proteomics, PTM, methylation, CNV, drug-response, calorimetry, and morphology extractors.
- Completed targeted and repo-wide validation sweeps.

## Validation status

- Targeted provenance-focused suite: `158 passed`
- Full repository pytest sweep: `303 passed`
