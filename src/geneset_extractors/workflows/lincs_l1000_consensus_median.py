"""Build LINCS L1000 chemical-perturbation median consensus signatures from GCTX."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from geneset_extractors.workflows.gtex_runtime_common import write_workflow_provenance_graph
from geneset_extractors.workflows.lincs_l1000_cp import MATRIX_PATH, ROW_ID_PATH, _h5py, _strings


PUBLIC_GCTX_URL = "https://lincs-dcic.s3.amazonaws.com/LINCS-sigs-2021/gctx/cd-coefficient/cp_coeff_mat.gctx"
SIGCOM_LIBRARY_UUID = "54198d6e-fe17-5ef8-91ac-02b425761653"
SIGCOM_METADATA_API = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
PERT_NAME_PATH = "0/META/COL/pert_name"


def inspect_gctx(gctx_path: Path) -> tuple[list[str], list[str], tuple[int, int]]:
    h5py = _h5py()
    with h5py.File(gctx_path, "r") as handle:
        missing = [path for path in (MATRIX_PATH, ROW_ID_PATH, PERT_NAME_PATH) if path not in handle]
        if missing:
            raise ValueError(f"GCTX is missing required dataset(s): {', '.join(missing)}")
        shape = tuple(int(value) for value in handle[MATRIX_PATH].shape)
        genes = [value.strip() for value in _strings(handle[ROW_ID_PATH][...])]
        pert_names = [value.strip() for value in _strings(handle[PERT_NAME_PATH][...])]
    if len(shape) != 2 or shape != (len(pert_names), len(genes)):
        raise ValueError(f"GCTX matrix shape {shape} is incompatible with {len(pert_names)} perturbation records and {len(genes)} row IDs")
    if not genes or any(not value for value in genes) or len(genes) != len(set(genes)):
        raise ValueError(f"{ROW_ID_PATH} must contain unique non-empty gene symbols")
    if not pert_names or any(not value for value in pert_names):
        raise ValueError(f"{PERT_NAME_PATH} must contain non-empty perturbagen names")
    return genes, pert_names, shape


def group_signature_indices(pert_names: list[str]) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = {}
    for index, pert_name in enumerate(pert_names):
        groups.setdefault(pert_name, []).append(index)
    return groups


def eligible_perturbagens(groups: dict[str, list[int]], min_signatures: int) -> list[str]:
    if min_signatures <= 0:
        raise ValueError("min_signatures must be positive")
    return sorted(name for name, indices in groups.items() if len(indices) >= min_signatures)


def partition_items(items: list[str], partition_index: int, partition_count: int) -> list[str]:
    if partition_count <= 0 or partition_index < 0 or partition_index >= partition_count:
        raise ValueError("partition_index must be in [0, partition_count)")
    start = len(items) * partition_index // partition_count
    end = len(items) * (partition_index + 1) // partition_count
    return items[start:end]


def rank_consensus(genes: list[str], coefficients: np.ndarray, top_n: int) -> tuple[list[str], list[str]]:
    if top_n <= 0 or len(genes) < top_n * 2:
        raise ValueError(f"Need at least {top_n * 2} genes; found {len(genes)}")
    values = np.asarray(coefficients, dtype=np.float64)
    if values.ndim != 1 or len(values) != len(genes) or not np.isfinite(values).all():
        raise ValueError("Invalid median CD-coefficient vector")
    order = np.lexsort((np.asarray(genes, dtype=str), -values))
    return [genes[index] for index in order[:top_n]], [genes[index] for index in order[-top_n:]]


def run(args) -> dict[str, object]:
    gctx_path = Path(args.gctx_path).resolve()
    if not gctx_path.is_file():
        raise FileNotFoundError(f"Missing cp_coeff_mat.gctx: {gctx_path}")
    top_n, min_signatures = int(args.top_n), int(args.min_signatures)
    genes, pert_names, shape = inspect_gctx(gctx_path)
    groups = group_signature_indices(pert_names)
    eligible = eligible_perturbagens(groups, min_signatures)
    partition_index, partition_count = int(args.partition_index), int(args.partition_count)
    selected = partition_items(eligible, partition_index, partition_count)
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    gmt_path = out_dir / "lincs_l1000_consensus_median.gmt"
    h5py = _h5py()
    with gmt_path.open("w", encoding="utf-8", newline="\n") as output, h5py.File(gctx_path, "r") as handle:
        matrix = handle[MATRIX_PATH]
        for pert_name in selected:
            values = np.asarray(matrix[groups[pert_name], :], dtype=np.float64)
            median = np.median(values, axis=0)
            up, down = rank_consensus(genes, median, top_n)
            output.write("\t".join([f"{pert_name} up", "", *up]) + "\n")
            output.write("\t".join([f"{pert_name} down", "", *down]) + "\n")
    manifest_path = out_dir / "lincs_l1000_consensus_median_partition.json"
    manifest = {"gctx_path": str(gctx_path), "public_gctx_url": PUBLIC_GCTX_URL, "sigcom_library_uuid": SIGCOM_LIBRARY_UUID, "sigcom_metadata_api": SIGCOM_METADATA_API, "matrix_shape": shape, "min_signatures": min_signatures, "top_n": top_n, "n_source_signatures": len(pert_names), "n_perturbagens": len(groups), "n_eligible_perturbagens": len(eligible), "partition_index": partition_index, "partition_count": partition_count, "perturbagens": selected, "n_generated_sets": len(selected) * 2}
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_workflow_provenance_graph(workflow_name="lincs_l1000_consensus_median", module_name="geneset_extractors.workflows.lincs_l1000_consensus_median", output_dir=out_dir, focus_output_path=gmt_path, output_paths=[(gmt_path, "consensus_median_gmt"), (manifest_path, "partition_manifest")], input_paths=[(gctx_path, "lincs_cp_coeff_mat_gctx")], parameters={"public_gctx_url": PUBLIC_GCTX_URL, "sigcom_library_uuid": SIGCOM_LIBRARY_UUID, "sigcom_metadata_api": SIGCOM_METADATA_API, "source_signature_schema": ["symbol", "CD-coefficient"], "grouping_key": "pert_name", "consensus_operation": "coordinate_wise_median", "min_signatures": min_signatures, "top_n": top_n, "ranking": "median CD-coefficient descending; gene symbol ascending", "partition_index": partition_index, "partition_count": partition_count, "n_source_signatures": len(pert_names), "n_eligible_perturbagens": len(eligible), "n_generated_sets": len(selected) * 2})
    return {"n_rows": len(selected) * top_n * 2, "n_sets": len(selected) * 2, "n_perturbagens": len(selected), "out_dir": str(out_dir)}


def merge_partitions(partition_dirs: list[Path], out_dir: Path) -> dict[str, int]:
    manifests = [json.loads((path / "lincs_l1000_consensus_median_partition.json").read_text(encoding="utf-8")) for path in partition_dirs]
    if not manifests:
        raise ValueError("At least one partition is required")
    partition_count = manifests[0]["partition_count"]
    if any(manifest["partition_count"] != partition_count for manifest in manifests):
        raise ValueError("Partition counts differ")
    if sorted(manifest["partition_index"] for manifest in manifests) != list(range(partition_count)):
        raise ValueError("Partitions must cover every index exactly once")
    out_dir.mkdir(parents=True, exist_ok=True)
    seen_terms: set[str] = set()
    with (out_dir / "lincs_l1000_consensus_median.gmt").open("w", encoding="utf-8", newline="\n") as output:
        for manifest, directory in sorted(zip(manifests, partition_dirs), key=lambda item: item[0]["partition_index"]):
            for line in (directory / "lincs_l1000_consensus_median.gmt").open(encoding="utf-8"):
                term = line.split("\t", 1)[0]
                if term in seen_terms:
                    raise ValueError(f"Duplicate output GMT term: {term}")
                seen_terms.add(term)
                output.write(line)
    return {"n_partitions": len(partition_dirs), "n_terms": len(seen_terms)}
