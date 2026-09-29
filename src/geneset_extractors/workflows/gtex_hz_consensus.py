"""GTEx V8 tissue, sex, and age consensus gene-set reconstruction.

This workflow owns the biological processing described in the transferred
legacy handoff.  The GTEx wrapper only selects HZ2 and invokes this module.
"""
from __future__ import annotations

import csv
import gzip
import math
import tempfile
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from geneset_extractors.core.provenance import activate_runtime_context
from geneset_extractors.core.metadata import write_metadata
from geneset_extractors.core.qc import write_run_summary_files
from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.workflows.gtex_runtime_common import (
    build_sample_metadata_rows,
    parse_gct_header,
    read_tsv,
    write_workflow_provenance_graph,
)


AGE_BINS = ("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
_NAME_TOKEN_RE = re.compile(r"[^A-Za-z0-9]+")


def _consensus_geneset_name(tissue: str, sex: str, age: str) -> str:
    tissue_token = _NAME_TOKEN_RE.sub("_", tissue).strip("_")
    return f"GTEx_Tissues_V8_Consensus_{tissue_token}_{sex}_{age}_up"


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open("r", encoding="utf-8", newline="")


def _sex_label(value: str) -> str | None:
    value = str(value).strip()
    if value in {"1", "1.0", "Male", "male", "MALE", "M"}:
        return "Male"
    if value in {"2", "2.0", "Female", "female", "FEMALE", "F"}:
        return "Female"
    return None


def _quantile_normalize_by_sample(values: np.ndarray) -> np.ndarray:
    """Quantile normalize columns, with stable ordering for tied expression."""
    order = np.argsort(values, axis=0, kind="mergesort")
    sorted_values = np.take_along_axis(values, order, axis=0)
    mean_by_rank = sorted_values.mean(axis=1)
    normalized = np.empty_like(values, dtype=float)
    for column in range(values.shape[1]):
        normalized[order[:, column], column] = mean_by_rank
    return normalized


def _sample_up_sets(gene_symbols: list[str], values: np.ndarray, cutoff: float) -> dict[str, set[str]]:
    """Return Harmonizome-style Up calls after robust per-sample standardization."""
    normalized = _quantile_normalize_by_sample(np.log2(np.maximum(values, 0.0) + 1.0))
    median = np.median(normalized, axis=0)
    mad = np.median(np.abs(normalized - median), axis=0)
    fallback = np.mean(np.abs(normalized - median), axis=0)
    scale = np.where(mad > 0, mad, fallback)
    scale = np.where(scale > 0, scale, 1.0)
    robust_z = (normalized - median) / scale
    up: dict[str, set[str]] = {}
    for column in range(robust_z.shape[1]):
        order = np.argsort(robust_z[:, column], kind="mergesort")
        ranks = np.empty(robust_z.shape[0], dtype=float)
        ranks[order] = np.arange(1, robust_z.shape[0] + 1, dtype=float)
        ecdf = ranks / float(robust_z.shape[0])
        up[str(column)] = {gene_symbols[index] for index in np.flatnonzero(ecdf >= cutoff)}
    return up


def _gct_row_count(expression_gct: Path) -> int:
    with _open_text(expression_gct) as handle:
        handle.readline()
        dimensions = handle.readline().strip().split("\t")
    try:
        return int(dimensions[0])
    except (IndexError, ValueError) as exc:
        raise ValueError("Expected a GCT dimensions line with a row count") from exc


def _write_sample_major_expression(
    expression_gct: Path,
    selected_ids: set[str],
    temp_dir: Path,
) -> tuple[list[str], list[str], np.memmap]:
    """Materialize only a disk-backed float32 expression matrix.

    The V8 TPM GCT is too large for nested Python lists or multiple dense
    in-memory arrays.  Gene-major layout permits efficient GCT ingestion; a
    sample-major copy makes the later rank computation contiguous per sample.
    """
    _, all_ids = parse_gct_header(expression_gct)
    sample_ids = [sample_id for sample_id in all_ids if sample_id in selected_ids]
    if not sample_ids:
        raise ValueError("No metadata-selected samples occur in the expression GCT")
    sample_positions = {sample_id: index + 2 for index, sample_id in enumerate(all_ids)}
    max_rows = _gct_row_count(expression_gct)
    gene_major = np.memmap(temp_dir / "gene_major.float32.mmap", mode="w+", dtype=np.float32, shape=(max_rows, len(sample_ids)))
    best_row_by_symbol: dict[str, tuple[int, float]] = {}
    retained_rows = 0
    with _open_text(expression_gct) as handle:
        handle.readline(); handle.readline()
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) < 3:
                continue
            symbol = str(row[1]).strip()
            if not symbol or symbol == "-":
                continue
            values = np.zeros(len(sample_ids), dtype=np.float32)
            for index, sample_id in enumerate(sample_ids):
                try:
                    values[index] = max(float(row[sample_positions[sample_id]] or 0.0), 0.0)
                except (IndexError, ValueError):
                    pass
            gene_major[retained_rows, :] = values
            value_sum = float(values.sum(dtype=np.float64))
            previous = best_row_by_symbol.get(symbol)
            # Preserve the established deterministic duplicate policy.
            if previous is None or value_sum > previous[1]:
                best_row_by_symbol[symbol] = (retained_rows, value_sum)
            retained_rows += 1
    if not best_row_by_symbol:
        raise ValueError("No usable gene-symbol rows were found in the expression GCT")
    symbols = sorted(best_row_by_symbol)
    source_rows = np.asarray([best_row_by_symbol[symbol][0] for symbol in symbols], dtype=np.intp)
    sample_major = np.memmap(temp_dir / "sample_major.float32.mmap", mode="w+", dtype=np.float32, shape=(len(sample_ids), len(symbols)))
    # Bounded blocks avoid a whole-matrix transpose allocation.
    for start in range(0, len(sample_ids), 128):
        stop = min(start + 128, len(sample_ids))
        sample_major[start:stop, :] = gene_major[source_rows, start:stop].T
    sample_major.flush()
    del gene_major
    return sample_ids, symbols, sample_major


def _ecdf_up_gene_indices(sample_values: np.ndarray, cutoff: float) -> np.ndarray:
    """Return the positive ECDF extreme without materializing transformed matrices.

    log2(TPM + 1), quantile normalization, and the positive-scale robust
    z-score are order-preserving within a sample.  The documented final Up call
    is solely ECDF rank >= cutoff, so this stable rank calculation is exactly
    equivalent to applying those intermediate monotonic transforms first.
    """
    order = np.argsort(sample_values, kind="mergesort")
    first_rank = int(math.ceil(float(cutoff) * len(order))) - 1
    return order[max(first_rank, 0):]


def _groups(sample_rows: list[dict[str, str]], expression_ids: set[str], min_samples: int) -> dict[tuple[str, str, str], list[str]]:
    grouped: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    for row in sample_rows:
        sample_id = row["sample_id"]
        sex = _sex_label(row.get("SEX", ""))
        age = row.get("age_bin", "")
        tissue = row.get("detailed_tissue", "").strip()
        if sample_id in expression_ids and sex and age in AGE_BINS and tissue:
            grouped[(tissue, sex, age)].append(sample_id)
    return {key: sorted(value) for key, value in grouped.items() if len(value) >= min_samples}


def run(args) -> dict[str, object]:
    expression_gct = Path(args.expression_gct).resolve()
    sample_attributes = Path(args.sample_attributes_tsv).resolve()
    subject_phenotypes = Path(args.subject_phenotypes_tsv).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    activate_runtime_context("gtex_hz_consensus", getattr(args, "provenance_overlay_json", None))

    _, gct_ids = parse_gct_header(expression_gct)
    metadata = build_sample_metadata_rows(
        sample_rows=read_tsv(sample_attributes), subject_rows=read_tsv(subject_phenotypes),
        gct_sample_ids=gct_ids, tissue_label="all detailed tissues", tissue_id="all_detailed_tissues",
        tissue_column=None, tissue_value=None, age_bins=list(AGE_BINS),
    )
    groups = _groups(metadata, set(gct_ids), int(args.min_samples_per_group))
    expected = getattr(args, "expected_group_count", None)
    if expected is not None and len(groups) != int(expected):
        raise ValueError(f"Expected {expected} tissue-sex-age groups; reconstructed {len(groups)}")
    gmt_path = out_dir / "genesets.gmt"
    support_path = out_dir / "gene_support.tsv"
    selected_path = out_dir / "geneset.tsv"
    full_path = out_dir / "geneset.full.tsv"
    ordered_groups = sorted(groups.items(), key=lambda item: f"{item[0][0]} {item[0][1]} {item[0][2]} Up")
    group_for_sample = {sample_id: group_index for group_index, (_, members) in enumerate(ordered_groups) for sample_id in members}
    selected_ids = set(group_for_sample)
    with tempfile.TemporaryDirectory(prefix="gtex_hz_consensus_", dir=out_dir) as temp_name:
        sample_ids, symbols, expression = _write_sample_major_expression(expression_gct, selected_ids, Path(temp_name))
        support_counts = np.zeros((len(ordered_groups), len(symbols)), dtype=np.uint16)
        for sample_index, sample_id in enumerate(sample_ids):
            up_indices = _ecdf_up_gene_indices(expression[sample_index, :], float(args.up_cutoff))
            support_counts[group_for_sample[sample_id], up_indices] += 1
        del expression
        gmt_sets: list[tuple[str, list[str]]] = []
        with support_path.open("w", encoding="utf-8", newline="") as support_handle, selected_path.open("w", encoding="utf-8", newline="") as selected_handle, full_path.open("w", encoding="utf-8", newline="") as full_handle:
            support_writer = csv.writer(support_handle, delimiter="\t", lineterminator="\n")
            support_writer.writerow(["term", "gene_symbol", "support_count", "support_fraction"])
            selected_writer = csv.writer(selected_handle, delimiter="\t", lineterminator="\n")
            full_writer = csv.writer(full_handle, delimiter="\t", lineterminator="\n")
            artifact_header = ["geneset_name", "gene_symbol", "support_count", "support_fraction"]
            selected_writer.writerow(artifact_header)
            full_writer.writerow(artifact_header)
            for group_index, ((tissue, sex, age), members) in enumerate(ordered_groups):
                counts = support_counts[group_index, :]
                eligible = np.flatnonzero(counts.astype(np.float64) / len(members) >= float(args.support_fraction))
                ranked = sorted(eligible, key=lambda index: (-(int(counts[index]) / len(members)), -int(counts[index]), symbols[index]))
                term = f"{tissue} {sex} {age} Up"
                geneset_name = _consensus_geneset_name(tissue, sex, age)
                genes = [symbols[index] for index in ranked[:int(args.top_n)]]
                gmt_sets.append((geneset_name, genes))
                for index in ranked:
                    row = [geneset_name, symbols[index], int(counts[index]), int(counts[index]) / len(members)]
                    support_writer.writerow([term, *row[1:]])
                    full_writer.writerow(row)
                for index in ranked[:int(args.top_n)]:
                    selected_writer.writerow([geneset_name, symbols[index], int(counts[index]), int(counts[index]) / len(members)])
        write_gmt(gmt_sets, gmt_path)
    graph = write_workflow_provenance_graph(
        workflow_name="gtex_hz_consensus", module_name=__name__, output_dir=out_dir,
        focus_output_path=gmt_path, output_paths=[(gmt_path, "gmt"), (support_path, "gene_support")],
        input_paths=[(expression_gct, "gtex_v8_tpm"), (sample_attributes, "sample_attributes_tsv_v8"), (subject_phenotypes, "subject_phenotypes_tsv_v8")],
        parameters={"grouping": "SMTSD x SEX x AGE; n_samples >= %d" % int(args.min_samples_per_group), "sample_signature": "log2(TPM+1), sample quantile normalization, robust median/MAD z-score with mean-absolute-deviation fallback, ECDF Up >= %.2f" % float(args.up_cutoff), "consensus": "support_fraction >= %.2f; top_n=%d" % (float(args.support_fraction), int(args.top_n)), "historical_drc_aggregation": "unavailable; scientifically comparable reconstruction, not set-equivalent"},
    )
    metadata_path = out_dir / "geneset.meta.json"
    write_metadata(metadata_path, {
        "schema_version": "1",
        "geneset_id": "gtex_hz2_consensus",
        "gene_set": {
            "id": "geneset:gtex_hz2_consensus",
            "name": "GTEx V8 tissue-sex-age consensus",
            "description": str(args.description),
        },
        "converter": {
            "name": "gtex_hz_consensus",
            "parameters": {"support_fraction": float(args.support_fraction), "top_n": int(args.top_n), "up_cutoff": float(args.up_cutoff)},
            "code": {"module": __name__},
            "execution": {"entrypoint": "geneset-extractors workflows gtex_hz_consensus"},
        },
        "input": {"files": [
            {"path": str(expression_gct), "role": "gtex_v8_tpm"},
            {"path": str(sample_attributes), "role": "sample_attributes_tsv_v8"},
            {"path": str(subject_phenotypes), "role": "subject_phenotypes_tsv_v8"},
        ]},
        "output": {"files": [
            {"path": "genesets.gmt", "role": "gmt"},
            {"path": "geneset.tsv", "role": "selected_program"},
            {"path": "geneset.full.tsv", "role": "full_scores"},
            {"path": "gene_support.tsv", "role": "gene_support"},
        ]},
        "provenance": {"focus_node_id": "geneset:gtex_hz2_consensus"},
        "library": "GTEx",
        "model_id": "HZ2",
        "model_family": "hz_consensus",
        "method": "gtex_hz_consensus",
        "description": str(args.description),
        "output_files": ["genesets.gmt", "geneset.tsv", "geneset.full.tsv", "gene_support.tsv"],
        "n_gene_sets": len(ordered_groups),
        "_provenance_overlay_json": getattr(args, "provenance_overlay_json", None),
        "_upstream_provenance_graph_path": str(graph),
    })
    (out_dir / "geneset.model.json").write_text(
        '{\n  "library": "GTEx",\n  "model_id": "HZ2",\n  "model_family": "hz_consensus",\n  "signature_pattern": "GTEx_Tissues_V8_Consensus_<tissue>_<sex>_<age>_up"\n}\n',
        encoding="utf-8",
    )
    write_run_summary_files(out_dir, {
        "workflow": "gtex_hz_consensus",
        "model_id": "HZ2",
        "n_groups": len(ordered_groups),
        "n_samples": len(sample_ids),
        "n_genes": len(symbols),
        "support_fraction": float(args.support_fraction),
        "top_n": int(args.top_n),
    })
    return {"out_dir": str(out_dir), "n_groups": len(groups), "n_samples": len(sample_ids)}
