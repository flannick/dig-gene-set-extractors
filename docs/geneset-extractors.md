# Gene Set Extractors Documentation Index

This is the repository-level documentation entrypoint.

For historical compatibility, many commands and modules still use the
`omics2geneset` name. New integrations can use the neutral
`geneset_extractors` namespace.

## Namespaces and CLI names

- Neutral namespace: `geneset_extractors` (new)
- Legacy namespace: `omics2geneset` (fully supported)
- CLI aliases:
  - `geneset-extractors`
  - `geneset_extractors`
  - `omics2geneset`

All three CLI names currently dispatch to the same implementation.

## Practical guides

- ATAC-seq: `docs/assays/atac/guide.md`
- RNA-seq: `docs/assays/rnaseq/guide.md`
- scRNA cNMF prep workflow: `docs/assays/rnaseq/scrna_cnmf_workflow.md`
- DNA methylation: `docs/assays/methylation/guide.md`
- DNA methylation resources: `docs/assays/methylation/resources.md`
- CNV segments: `docs/assays/cnv/guide.md`

## Methods notes

- ATAC-seq methods: `docs/assays/atac/methods.tex`
- RNA-seq methods: `docs/assays/rnaseq/methods.tex`
- DNA methylation methods: `docs/assays/methylation/methods.tex`
- CNV methods: `docs/assays/cnv/methods.tex`
- Methods index: `docs/methods.tex`

Compatibility aliases remain at historical paths such as `docs/atac-seq2geneset.md` and `docs/rna-seq_methods.tex`.
