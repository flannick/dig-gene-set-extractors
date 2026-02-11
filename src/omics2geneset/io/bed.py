from __future__ import annotations

from pathlib import Path


def read_bed(path: str | Path) -> list[dict[str, object]]:
    peaks: list[dict[str, object]] = []
    with Path(path).open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            peaks.append(
                {
                    "chrom": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "columns": parts,
                }
            )
    return peaks


def read_peak_weights_tsv(path: str | Path) -> dict[tuple[str, int, int], float]:
    out: dict[tuple[str, int, int], float] = {}
    with Path(path).open("r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        for key in ("chrom", "start", "end", "weight"):
            if key not in index:
                raise ValueError(f"peak_weights_tsv missing required column: {key}")
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            out[(p[index["chrom"]], int(p[index["start"]]), int(p[index["end"]]))] = float(p[index["weight"]])
    return out
