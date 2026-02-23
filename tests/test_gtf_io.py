import gzip
from pathlib import Path

from omics2geneset.io.gtf import read_genes_from_gtf


def test_read_genes_from_gtf_accepts_gzip(tmp_path: Path):
    gtf_gz = tmp_path / "toy.gtf.gz"
    with gzip.open(gtf_gz, "wt", encoding="utf-8") as fh:
        fh.write(
            "\n".join(
                [
                    "##gtf-version 2.2",
                    'chr1\tsrc\tgene\t101\t200\t.\t+\t.\tgene_id "ENSG000001"; gene_name "GENE1"; gene_type "protein_coding";',
                    'chr2\tsrc\tgene\t301\t500\t.\t-\t.\tgene_id "ENSG000002"; gene_name "GENE2"; gene_type "protein_coding";',
                    "",
                ]
            )
        )

    genes = read_genes_from_gtf(gtf_gz)
    assert [g.gene_id for g in genes] == ["ENSG000001", "ENSG000002"]
    assert genes[0].chrom == "chr1"
    assert genes[0].tss == 100
    assert genes[1].chrom == "chr2"
    assert genes[1].tss == 499
