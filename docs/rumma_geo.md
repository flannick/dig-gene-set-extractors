# RummaGEO reconstruction

## Normal workflow

Run both user-facing libraries with one command:

```bash
geneset-extractors convert rumma_geo_all --human_gmt human-geo-auto.gmt.gz --mouse_gmt mouse-geo-auto.gmt.gz --human_gene_info Homo_sapiens.gene_info.gz --mouse_gene_info Mus_musculus.gene_info.gz --gene_orthologs gene_orthologs.gz --out_dir output
```

It caches SigCom LINCS metadata, derives deterministic `pert_name` drug terms,
acquires and caches RummaGEO GraphQL records, creates selection manifests, and
writes `output/genesets/all_signatures/models/HZ2/extractor/` and
`output/genesets/all_signatures/models/HZ1/extractor/`; each model's staged
artifacts live under its sibling `workflow/` directory. Reuse caches by default; `--refresh_sources` explicitly
reacquires mutable upstream resources. The SigCom URL is
`https://s3.dev.maayanlab.cloud/sigcom-lincs/ranker/signatures_meta.json`, as
specified by `RummaGEODrug.ipynb` at HarmonizomePythonScripts commit
`965d3a7299cdeaa8d54740b31093b80cebd5523b`. Its cached bytes and SHA-256 are
the reproducible input; the URL is not an immutable historical snapshot.

Supply the two GMTs plus pinned human/mouse NCBI gene-info and ortholog files.
Current RummaGEO, SigCom, and NCBI resources can differ from unavailable
historical snapshots, so this is a method-faithful reproducible reconstruction,
not a claim of byte-for-byte historical reproduction. Recovered snapshots can
be substituted without changing the algorithm.

## Advanced workflow

For offline inspection or staged reruns, use `rumma_geo_acquire`, then
`rumma_geo_selection`, then `rumma_geo` directly. `rumma_geo_acquire` records
GraphQL pagination and cached-query checksums; selection records control and
reversed classification; reconstruction consumes the resulting
`selection_manifest.tsv`.

`rumma_geo` reconstructs one of two RummaGEO models: `HZ1` (drug perturbations)
or `HZ2` (gene perturbations). It ports the downstream behavior of the Ma'ayan Lab
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
