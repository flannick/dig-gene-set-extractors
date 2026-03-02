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

- ATAC-seq: `docs/atac-seq2geneset.md`
- RNA-seq: `docs/rna-seq-to-geneset.md`
- scRNA cNMF prep workflow: `docs/scrna_cnmf_workflow.md`
- DNA methylation: `docs/methylation2geneset.md`
- CNV segments: `docs/cnv2geneset.md`

## Methods notes

- ATAC-seq methods: `docs/atac-seq_methods.tex`
- RNA-seq methods: `docs/rna-seq_methods.tex`
- DNA methylation methods: `docs/methylation_methods.tex`
- CNV methods: `docs/cnv_methods.tex`
- Methods index: `docs/methods.tex`
