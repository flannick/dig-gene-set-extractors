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

    @staticmethod
    def _row_score(row: dict[str, object]) -> float:
        if "score" in row and row["score"] is not None:
            return float(row["score"])
        if "weight" in row and row["weight"] is not None:
            return float(row["weight"])
        return 0.0

    def sort_desc(self) -> "GeneWeights":
        self.rows = sorted(self.rows, key=lambda r: (-self._row_score(r), str(r["gene_id"])))
        for idx, row in enumerate(self.rows, start=1):
            row["rank"] = idx
        return self

    def l1_normalize(self) -> "GeneWeights":
        total = sum(self._row_score(r) for r in self.rows)
        if total <= 0:
            return self
        for row in self.rows:
            v = self._row_score(row) / total
            if "score" in row:
                row["score"] = v
            if "weight" in row:
                row["weight"] = v
            if "score" not in row and "weight" not in row:
                row["score"] = v
        return self

    def to_tsv(self, path: str | Path) -> None:
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = ["gene_id", "score"]
        if any("weight" in r for r in self.rows):
            fieldnames.append("weight")
        if any("gene_symbol" in r and r["gene_symbol"] is not None for r in self.rows):
            fieldnames.append("gene_symbol")
        if any("rank" in r for r in self.rows):
            fieldnames.append("rank")
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for row in self.rows:
                out_row = dict(row)
                if "score" not in out_row:
                    out_row["score"] = self._row_score(row)
                writer.writerow(out_row)
