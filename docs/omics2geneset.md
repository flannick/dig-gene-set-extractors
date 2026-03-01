# omics2geneset Documentation Index

`omics2geneset` now has assay-specific practical guides and methods notes.

Use this page as the entrypoint, then jump to the assay you are running.

## Practical guides

- ATAC-seq: `docs/atac-seq2geneset.md`
- RNA-seq: `docs/rna-seq-to-geneset.md`
- DNA methylation: `docs/methylation2geneset.md`

## Methods notes

- ATAC-seq methods: `docs/atac-seq_methods.tex`
- RNA-seq methods: `docs/rna-seq_methods.tex`
- DNA methylation methods: `docs/methylation_methods.tex`
- Methods index: `docs/methods.tex`

## Choosing a guide

- If your input is peaks, peak matrices, or scATAC 10x matrices: start with `docs/atac-seq2geneset.md`.
- If your input is RNA differential expression results or external scRNA program loadings (for example cNMF/scHPF gene loadings): start with `docs/rna-seq-to-geneset.md`.
- If your input is differential CpG/probe or DMR tables: start with `docs/methylation2geneset.md`.

## Backward compatibility note

Older references to `docs/omics2geneset.md` and `docs/methods.tex` remain valid.
Those files now act as index pages pointing to assay-specific documentation.
