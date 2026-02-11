from __future__ import annotations

from collections import defaultdict
import gzip
from pathlib import Path
import re


FEATURE_COORD_RE = re.compile(r"^([^:]+):(\d+)-(\d+)$")


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _resolve_one(base: Path, names: list[str]) -> Path:
    for name in names:
        p = base / name
        if p.exists():
            return p
    raise FileNotFoundError(f"None of expected files found in {base}: {names}")


def _parse_feature_line(parts: list[str]) -> dict[str, object]:
    if len(parts) >= 3 and parts[1].isdigit() and parts[2].isdigit():
        return {"chrom": parts[0], "start": int(parts[1]), "end": int(parts[2])}

    if parts:
        m = FEATURE_COORD_RE.match(parts[0])
        if m:
            chrom, start, end = m.groups()
            return {"chrom": chrom, "start": int(start), "end": int(end)}

    raise ValueError(f"Could not parse peak coordinates from feature row: {'\t'.join(parts)}")


def read_10x_matrix_dir(
    matrix_dir: str | Path,
) -> tuple[list[dict[str, object]], list[str], list[tuple[int, int, float]], dict[str, Path]]:
    base = Path(matrix_dir)
    matrix_path = _resolve_one(base, ["matrix.mtx", "matrix.mtx.gz"])
    barcodes_path = _resolve_one(base, ["barcodes.tsv", "barcodes.tsv.gz"])
    peaks_path = _resolve_one(base, ["peaks.bed", "features.tsv", "features.tsv.gz"])

    peaks: list[dict[str, object]] = []
    with _open_text(peaks_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            peaks.append(_parse_feature_line(p))

    barcodes: list[str] = []
    with _open_text(barcodes_path) as fh:
        for line in fh:
            if line.strip():
                barcodes.append(line.strip().split("\t")[0])

    entries: list[tuple[int, int, float]] = []
    with _open_text(matrix_path) as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            dims = line.strip().split()
            if len(dims) == 3:
                break
        for line in fh:
            if not line.strip():
                continue
            r, c, v = line.strip().split()
            entries.append((int(r) - 1, int(c) - 1, float(v)))

    files = {
        "matrix": matrix_path,
        "barcodes": barcodes_path,
        "peaks_or_features": peaks_path,
    }
    return peaks, barcodes, entries, files


def summarize_peaks(entries: list[tuple[int, int, float]], n_peaks: int, n_cells: int, cell_indices: list[int], method: str) -> list[float]:
    sums = [0.0] * n_peaks
    nonzero = [0] * n_peaks
    selected = set(cell_indices)
    for peak_i, cell_i, val in entries:
        if cell_i not in selected:
            continue
        sums[peak_i] += val
        if val != 0:
            nonzero[peak_i] += 1
    denom = max(len(cell_indices), 1)
    if method == "sum_counts":
        return sums
    if method == "mean_counts":
        return [x / denom for x in sums]
    if method == "frac_cells_nonzero":
        return [x / denom for x in nonzero]
    raise ValueError(f"Unsupported peak_summary: {method}")


def read_groups_tsv(path: str | Path) -> dict[str, str]:
    groups: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        if "barcode" not in idx or "group" not in idx:
            raise ValueError("groups_tsv must contain barcode and group columns")
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            groups[p[idx["barcode"]]] = p[idx["group"]]
    return groups


def make_group_indices(barcodes: list[str], groups: dict[str, str] | None) -> dict[str, list[int]]:
    if not groups:
        return {"all": list(range(len(barcodes)))}
    by_group: dict[str, list[int]] = defaultdict(list)
    for i, bc in enumerate(barcodes):
        g = groups.get(bc)
        if g is not None:
            by_group[g].append(i)
    return dict(by_group)
