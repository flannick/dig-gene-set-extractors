from __future__ import annotations

import csv
from pathlib import Path


def read_region_gene_links_tsv(path: str | Path) -> list[dict[str, object]]:
    links: list[dict[str, object]] = []
    with Path(path).open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"chrom", "start", "end", "gene_id", "link_weight"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(
                "region_gene_links_tsv missing required columns: " + ", ".join(missing)
            )
        for row in reader:
            chrom = str(row["chrom"]).strip()
            gene_id = str(row["gene_id"]).strip()
            if not chrom or not gene_id:
                continue
            start = int(row["start"])
            end = int(row["end"])
            if end <= start:
                continue
            link_weight = float(row["link_weight"])
            links.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "gene_id": gene_id,
                    "link_weight": link_weight,
                }
            )
    return links
