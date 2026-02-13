from __future__ import annotations

from collections import defaultdict
import gzip
from pathlib import Path
import re
from typing import Iterator


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

    joined = "\t".join(parts)
    raise ValueError(f"Could not parse peak coordinates from feature row: {joined}")


def read_10x_matrix_dir(
    matrix_dir: str | Path,
) -> tuple[list[dict[str, object]], list[str], dict[str, Path]]:
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

    files = {
        "matrix": matrix_path,
        "barcodes": barcodes_path,
        "peaks_or_features": peaks_path,
    }
    return peaks, barcodes, files


def iter_mtx_entries(matrix_path: Path) -> Iterator[tuple[int, int, float]]:
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
            yield int(r) - 1, int(c) - 1, float(v)


def summarize_peaks(matrix_path: Path, n_peaks: int, cell_indices: list[int], method: str) -> list[float]:
    sums = [0.0] * n_peaks
    nonzero = [0] * n_peaks
    selected = set(cell_indices)
    for peak_i, cell_i, val in iter_mtx_entries(matrix_path):
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


def summarize_peaks_by_group(
    matrix_path: Path,
    n_peaks: int,
    group_indices: dict[str, list[int]],
    method: str,
) -> dict[str, list[float]]:
    group_names = list(group_indices)
    cell_to_group_idx: dict[int, int] = {}
    for gi, group in enumerate(group_names):
        for cell_idx in group_indices[group]:
            cell_to_group_idx[cell_idx] = gi

    sums = [[0.0] * n_peaks for _ in group_names]
    nonzero = [[0] * n_peaks for _ in group_names]

    for peak_i, cell_i, val in iter_mtx_entries(matrix_path):
        gi = cell_to_group_idx.get(cell_i)
        if gi is None:
            continue
        sums[gi][peak_i] += val
        if val != 0:
            nonzero[gi][peak_i] += 1

    result: dict[str, list[float]] = {}
    for gi, group in enumerate(group_names):
        denom = max(len(group_indices[group]), 1)
        if method == "sum_counts":
            result[group] = sums[gi]
        elif method == "mean_counts":
            result[group] = [x / denom for x in sums[gi]]
        elif method == "frac_cells_nonzero":
            result[group] = [x / denom for x in nonzero[gi]]
        else:
            raise ValueError(f"Unsupported peak_summary: {method}")
    return result


def summarize_group_vs_rest(
    matrix_path: Path,
    n_peaks: int,
    n_cells: int,
    group_indices: dict[str, list[int]],
    method: str,
) -> tuple[dict[str, list[float]], dict[str, list[float]], dict[str, int]]:
    group_names = list(group_indices)
    cell_to_group_idx: dict[int, int] = {}
    for gi, group in enumerate(group_names):
        for cell_idx in group_indices[group]:
            cell_to_group_idx[cell_idx] = gi

    group_sums = [[0.0] * n_peaks for _ in group_names]
    group_nonzero = [[0] * n_peaks for _ in group_names]
    all_sums = [0.0] * n_peaks
    all_nonzero = [0] * n_peaks

    for peak_i, cell_i, val in iter_mtx_entries(matrix_path):
        all_sums[peak_i] += val
        if val != 0:
            all_nonzero[peak_i] += 1
        gi = cell_to_group_idx.get(cell_i)
        if gi is None:
            continue
        group_sums[gi][peak_i] += val
        if val != 0:
            group_nonzero[gi][peak_i] += 1

    group_summary: dict[str, list[float]] = {}
    rest_summary: dict[str, list[float]] = {}
    group_sizes = {g: len(group_indices[g]) for g in group_names}

    for gi, group in enumerate(group_names):
        n_group = max(group_sizes[group], 1)
        n_rest = max(n_cells - group_sizes[group], 1)
        if method == "sum_counts":
            g_vals = group_sums[gi]
            r_vals = [all_sums[i] - group_sums[gi][i] for i in range(n_peaks)]
        elif method == "mean_counts":
            g_vals = [x / n_group for x in group_sums[gi]]
            r_vals = [(all_sums[i] - group_sums[gi][i]) / n_rest for i in range(n_peaks)]
        elif method == "frac_cells_nonzero":
            g_vals = [x / n_group for x in group_nonzero[gi]]
            r_vals = [(all_nonzero[i] - group_nonzero[gi][i]) / n_rest for i in range(n_peaks)]
        else:
            raise ValueError(f"Unsupported peak_summary: {method}")
        group_summary[group] = g_vals
        rest_summary[group] = r_vals
    return group_summary, rest_summary, group_sizes


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
