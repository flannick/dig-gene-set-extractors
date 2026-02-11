from __future__ import annotations

from dataclasses import dataclass
import csv
from pathlib import Path


@dataclass(frozen=True)
class GenomicInterval:
    chrom: str
    start: int
    end: int


@dataclass(frozen=True)
class Gene:
    gene_id: str
    gene_symbol: str | None
    chrom: str
    tss: int
    strand: str
    gene_start: int
    gene_end: int


class GeneWeights:
    def __init__(self, rows: list[dict[str, object]]):
        self.rows = rows

    def sort_desc(self) -> "GeneWeights":
        self.rows = sorted(self.rows, key=lambda r: float(r["weight"]), reverse=True)
        for idx, row in enumerate(self.rows, start=1):
            row["rank"] = idx
        return self

    def l1_normalize(self) -> "GeneWeights":
        total = sum(float(r["weight"]) for r in self.rows)
        if total <= 0:
            return self
        for row in self.rows:
            row["weight"] = float(row["weight"]) / total
        return self

    def to_tsv(self, path: str | Path) -> None:
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = ["gene_id", "weight"]
        if any("gene_symbol" in r and r["gene_symbol"] is not None for r in self.rows):
            fieldnames.append("gene_symbol")
        if any("rank" in r for r in self.rows):
            fieldnames.append("rank")
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for row in self.rows:
                writer.writerow(row)
