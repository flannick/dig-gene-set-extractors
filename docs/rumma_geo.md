# RummaGEO reconstruction

`rumma_geo` reconstructs one of two RummaGEO models: `gene_perturbations` or
`drug_perturbations`. It ports the downstream behavior of the Ma'ayan Lab
notebooks: exact human protein-coding symbol filtering, mouse-symbol to human
ortholog conversion, complete removal of duplicate term/gene pairs, and a
five-gene minimum for directional GMT sets. It intentionally performs no
synonym rescue.

The notebook expression for a reversed comparison mapped both directions to
`dn`, an implementation bug. The standardized RummaGEO GMT correctly swaps
`up` and `dn`; this converter follows that production behavior.

The historical notebooks first queried Harmonizome for selection rows. Those
query responses are not derivable faithfully from a raw GEO GMT. Supply them
as `selection_manifest` with `source_term`, `model_id`, and `status`; use
`normalized_term` to preserve the exact notebook-derived label. Rows with
`status` other than `signature` or `reversed` are excluded. For drug models,
the historical false-positive search terms are excluded.

The required source manifest is JSON with one URL and version for each input:

```json
{"sources":{"human_rummageo_gmt":{"url":"https://…","version":"2024-11"},"mouse_rummageo_gmt":{"url":"https://…","version":"2024-11"},"recorded_selection_manifest":{"url":"https://…","version":"2024-11"},"ncbi_human_gene_info":{"url":"https://…","version":"2024-11-01"},"ncbi_mouse_gene_info":{"url":"https://…","version":"2024-11-01"},"ncbi_gene_orthologs":{"url":"https://…","version":"2024-11-01"}}}
```

The emitted metadata records these pinned identifiers and local SHA-256 hashes.
Use `--legacy_gmt` only to calculate validation counts, precision, recall,
Jaccard, and per-set Jaccard in `reconstruction_diagnostics.json`; it cannot
alter reconstruction output.
