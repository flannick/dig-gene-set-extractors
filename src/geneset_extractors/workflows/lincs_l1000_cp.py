"""Stream LINCS L1000 chemical-perturbation CD signatures from a GCTX file."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable

import numpy as np

from geneset_extractors.workflows.gtex_runtime_common import write_workflow_provenance_graph

PUBLIC_GCTX_URL = "https://lincs-dcic.s3.amazonaws.com/LINCS-sigs-2021/gctx/cd-coefficient/cp_coeff_mat.gctx"
MATRIX_PATH = "0/DATA/0/matrix"
ROW_ID_PATH = "0/META/ROW/id"
LINCS_ID_PATH = "0/META/COL/lincs_id"
CELL_LINE_PATH = "0/META/COL/cell_line"
PERT_TIME_PATH = "0/META/COL/pert_time"


def _h5py():
    try:
        import h5py
    except ImportError as error:
        raise RuntimeError("lincs_l1000_cp requires h5py to read GCTX inputs") from error
    return h5py


def _strings(values: Iterable[object]) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def inspect_gctx(gctx_path: Path) -> tuple[list[str], list[str], tuple[int, int]]:
    """Validate the required public GCTX datasets and their dimensions."""
    h5py = _h5py()
    with h5py.File(gctx_path, "r") as handle:
        missing = [path for path in (MATRIX_PATH, ROW_ID_PATH, LINCS_ID_PATH) if path not in handle]
        if missing:
            raise ValueError(f"GCTX is missing required dataset(s): {', '.join(missing)}")
        shape = tuple(int(value) for value in handle[MATRIX_PATH].shape)
        genes = [value.strip() for value in _strings(handle[ROW_ID_PATH][...])]
        signatures = [value.strip() for value in _strings(handle[LINCS_ID_PATH][...])]
    if len(shape) != 2 or shape != (len(signatures), len(genes)):
        raise ValueError(f"GCTX matrix shape {shape} is incompatible with {len(signatures)} lincs_id values and {len(genes)} row IDs")
    if not genes or any(not value for value in genes) or len(genes) != len(set(genes)):
        raise ValueError(f"{ROW_ID_PATH} must contain unique non-empty gene symbols")
    if not signatures or any(not value for value in signatures):
        raise ValueError(f"{LINCS_ID_PATH} must contain non-empty lincs_id values")
    return genes, signatures, shape


def resolve_last_occurrences(signatures: list[str]) -> tuple[list[int], dict[str, list[int]]]:
    """Return retained raw-column indices and duplicate groups; later columns win."""
    positions: dict[str, list[int]] = {}
    for index, signature in enumerate(signatures):
        positions.setdefault(signature, []).append(index)
    retained = sorted(indices[-1] for indices in positions.values())
    duplicates = {signature: indices for signature, indices in positions.items() if len(indices) > 1}
    return retained, duplicates


def duplicate_vector_counts(matrix, duplicate_groups: dict[str, list[int]]) -> tuple[int, int]:
    exact, differing = 0, 0
    for indices in duplicate_groups.values():
        replacement = matrix[indices[-1], :]
        if all(np.array_equal(matrix[index, :], replacement) for index in indices[:-1]):
            exact += 1
        else:
            differing += 1
    return exact, differing


def resolve_range(n_signatures: int, start_index: int, end_index: int | None) -> tuple[int, int]:
    end = n_signatures if end_index is None else end_index
    if start_index < 0 or end < start_index or end > n_signatures:
        raise ValueError(f"Invalid signature range [{start_index}, {end}) for {n_signatures} signatures")
    return start_index, end


def plan_cell_time_partitions(gctx_path: Path, out_dir: Path, max_signatures_per_task: int) -> int:
    """Write deterministic retained-signature worklists grouped by cell line/time."""
    if max_signatures_per_task <= 0:
        raise ValueError("max_signatures_per_task must be positive")
    _genes, signatures, _shape = inspect_gctx(gctx_path)
    retained_indices, _duplicates = resolve_last_occurrences(signatures)
    h5py = _h5py()
    with h5py.File(gctx_path, "r") as handle:
        missing = [path for path in (CELL_LINE_PATH, PERT_TIME_PATH) if path not in handle]
        if missing:
            raise ValueError(f"GCTX is missing partition metadata dataset(s): {', '.join(missing)}")
        cell_lines = [value.strip() or "NA" for value in _strings(handle[CELL_LINE_PATH][...])]
        pert_times = [value.strip() or "NA" for value in _strings(handle[PERT_TIME_PATH][...])]
    if len(cell_lines) != len(signatures) or len(pert_times) != len(signatures):
        raise ValueError("GCTX partition metadata length is incompatible with lincs_id")
    groups: dict[tuple[str, str], list[int]] = {}
    for raw_index in retained_indices:
        groups.setdefault((cell_lines[raw_index], pert_times[raw_index]), []).append(raw_index)
    indices_dir = out_dir / "indices"
    indices_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    task_number = 0
    for cell_line, pert_time in sorted(groups):
        indices = groups[(cell_line, pert_time)]
        for offset in range(0, len(indices), max_signatures_per_task):
            task_number += 1
            selected = indices[offset : offset + max_signatures_per_task]
            task_id = f"hz4_{task_number:04d}"
            index_path = indices_dir / f"{task_id}.tsv"
            index_path.write_text("raw_index\n" + "".join(f"{index}\n" for index in selected), encoding="utf-8")
            rows.append({"task_id": task_id, "cell_line": cell_line, "pert_time": pert_time, "raw_indices_tsv": str(index_path), "n_signatures": len(selected)})
    with (out_dir / "task_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["task_id", "cell_line", "pert_time", "raw_indices_tsv", "n_signatures"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def _read_raw_indices(path_text: str, retained_indices: list[int]) -> list[int]:
    path = Path(path_text).resolve()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "raw_index" not in reader.fieldnames:
            raise ValueError(f"Index file lacks raw_index column: {path}")
        indices = [int(row["raw_index"]) for row in reader]
    if not indices or len(indices) != len(set(indices)) or set(indices).difference(retained_indices):
        raise ValueError(f"Index file must contain unique retained GCTX column indices: {path}")
    return indices


def rank_signature(genes: list[str], coefficients: np.ndarray, top_n: int) -> tuple[list[str], list[str]]:
    if len(genes) < top_n * 2:
        raise ValueError(f"Need at least {top_n * 2} genes; found {len(genes)}")
    values = np.asarray(coefficients, dtype=np.float64)
    if values.ndim != 1 or len(values) != len(genes) or not np.isfinite(values).all():
        raise ValueError("Invalid CD-coefficient vector")
    order = np.lexsort((np.asarray(genes, dtype=str), -values))
    return [genes[index] for index in order[:top_n]], [genes[index] for index in order[-top_n:]]


def _write_gmt(path: Path, rows: list[tuple[str, str, list[str]]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for term, direction, genes in rows:
            handle.write("\t".join([f"{term} {direction}", "", *genes]) + "\n")


def _write_signed(path: Path, rows: list[tuple[str, str, list[str]]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["term", "gene_id", "gene_symbol", "score", "sign"], lineterminator="\n")
        writer.writeheader()
        for term, direction, genes in rows:
            for rank, gene in enumerate(genes, start=1):
                writer.writerow({"term": term, "gene_id": gene, "gene_symbol": gene, "score": len(genes) - rank + 1, "sign": 1 if direction == "up" else -1})


def run(args) -> dict[str, object]:
    gctx_path = Path(args.gctx_path).resolve()
    if not gctx_path.is_file():
        raise FileNotFoundError(f"Missing cp_coeff_mat.gctx: {gctx_path}")
    top_n, block_size = int(args.top_n), int(args.block_size)
    if top_n <= 0 or block_size <= 0:
        raise ValueError("top_n and block_size must be positive")
    genes, signatures, shape = inspect_gctx(gctx_path)
    retained_indices, duplicate_groups = resolve_last_occurrences(signatures)
    raw_indices_tsv = getattr(args, "raw_indices_tsv", None)
    if raw_indices_tsv:
        selected_indices = _read_raw_indices(raw_indices_tsv, retained_indices)
        start, end = None, None
    else:
        start, end = resolve_range(len(retained_indices), int(args.start_index), args.end_index)
        selected_indices = retained_indices[start:end]
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    n_sets = 0
    h5py = _h5py()
    gmt_path, signed_path = out_dir / "l1000_cp.gmt", out_dir / "lincs_l1000_cp_signed_term_gene.tsv"
    with gmt_path.open("w", encoding="utf-8", newline="\n") as gmt_handle, signed_path.open("w", encoding="utf-8", newline="") as signed_handle, h5py.File(gctx_path, "r") as handle:
        writer = csv.DictWriter(signed_handle, delimiter="\t", fieldnames=["term", "gene_id", "gene_symbol", "score", "sign"], lineterminator="\n")
        writer.writeheader()
        matrix = handle[MATRIX_PATH]
        exact_duplicate_groups, differing_duplicate_groups = duplicate_vector_counts(matrix, duplicate_groups)
        for block_start in range(0, len(selected_indices), block_size):
            raw_indices = selected_indices[block_start : block_start + block_size]
            for raw_index in raw_indices:
                coefficients = matrix[raw_index, :]
                term = signatures[raw_index]
                up, down = rank_signature(genes, coefficients, top_n)
                for direction, selected in (("up", up), ("down", down)):
                    gmt_handle.write("\t".join([f"{term} {direction}", "", *selected]) + "\n")
                    for rank, gene in enumerate(selected, start=1):
                        writer.writerow({"term": term, "gene_id": gene, "gene_symbol": gene, "score": top_n - rank + 1, "sign": 1 if direction == "up" else -1})
                    n_sets += 1
    manifest_path = out_dir / "lincs_l1000_cp_partition.json"
    duplicate_summary = {"raw_signature_columns": len(signatures), "unique_lincs_id": len(retained_indices), "duplicate_lincs_id_groups": len(duplicate_groups), "exact_duplicate_vector_groups": exact_duplicate_groups, "differing_duplicate_vector_groups": differing_duplicate_groups, "resolution_policy": "last GCTX column occurrence wins"}
    manifest_path.write_text(json.dumps({"gctx_path": str(gctx_path), "public_url": PUBLIC_GCTX_URL, "matrix_shape": shape, "start_index": start, "end_index": end, "raw_indices_tsv": raw_indices_tsv, "n_signatures": len(selected_indices), "n_sets": n_sets, "top_n": top_n, "duplicate_resolution": duplicate_summary}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_workflow_provenance_graph(workflow_name="lincs_l1000_cp", module_name="geneset_extractors.workflows.lincs_l1000_cp", output_dir=out_dir, focus_output_path=signed_path, output_paths=[(signed_path, "signed_term_gene_tsv"), (gmt_path, "per_signature_gmt"), (manifest_path, "partition_manifest")], input_paths=[(gctx_path, "lincs_cp_coeff_mat_gctx")], parameters={"public_gctx_url": PUBLIC_GCTX_URL, "required_datasets": [MATRIX_PATH, ROW_ID_PATH, LINCS_ID_PATH], "top_n": top_n, "ranking": "CD-coefficient descending; symbol ascending", "start_index": start, "end_index": end, "raw_indices_tsv": raw_indices_tsv, "n_source_signatures": len(signatures), "n_unique_lincs_id": len(retained_indices), "duplicate_resolution": duplicate_summary, "n_generated_sets": n_sets})
    return {"n_rows": n_sets * top_n, "n_sets": n_sets, "out_dir": str(out_dir)}


def merge_partitions(partition_dirs: list[Path], out_dir: Path) -> dict[str, int]:
    manifests = [json.loads((path / "lincs_l1000_cp_partition.json").read_text(encoding="utf-8")) for path in partition_dirs]
    ordered = sorted(zip(manifests, partition_dirs), key=lambda item: item[0]["start_index"])
    expected_start, seen_terms = 0, set()
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "l1000_cp.gmt").open("w", encoding="utf-8", newline="\n") as output:
        for manifest, directory in ordered:
            if manifest["start_index"] != expected_start:
                raise ValueError("Partition ranges are not contiguous from signature index 0")
            expected_start = manifest["end_index"]
            for line in (directory / "l1000_cp.gmt").open(encoding="utf-8"):
                term = line.split("\t", 1)[0]
                if term in seen_terms:
                    raise ValueError(f"Duplicate output GMT term: {term}")
                seen_terms.add(term)
                output.write(line)
    return {"n_partitions": len(ordered), "n_terms": len(seen_terms), "end_index": expected_start}
