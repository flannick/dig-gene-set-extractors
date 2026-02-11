from __future__ import annotations

from pathlib import Path

from omics2geneset.core.models import Gene


def _parse_attrs(attr_field: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for chunk in attr_field.split(";"):
        c = chunk.strip()
        if not c:
            continue
        if " " not in c:
            continue
        k, v = c.split(" ", 1)
        attrs[k] = v.strip().strip('"')
    return attrs


def read_genes_from_gtf(path: str | Path, gene_id_field: str = "gene_id") -> list[Gene]:
    genes: list[Gene] = []
    with Path(path).open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, strand, _, attrs_str = line.rstrip("\n").split("\t")
            if feature != "gene":
                continue
            attrs = _parse_attrs(attrs_str)
            gene_id = attrs.get(gene_id_field)
            if not gene_id:
                continue
            gene_symbol = attrs.get("gene_name")
            start_1 = int(start)
            end_1 = int(end)
            start_0 = start_1 - 1
            end_0 = end_1
            tss = start_0 if strand == "+" else end_0 - 1
            genes.append(
                Gene(
                    gene_id=gene_id,
                    gene_symbol=gene_symbol,
                    chrom=chrom,
                    tss=tss,
                    strand=strand,
                    gene_start=start_0,
                    gene_end=end_0,
                )
            )
    return genes
