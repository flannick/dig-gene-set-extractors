from __future__ import annotations

import csv
import gzip
from pathlib import Path
import sys
from typing import Any, Iterable

from geneset_extractors.core.metadata import _resolve_git_commit, current_invocation_context
from geneset_extractors.core.provenance import (
    REPO_URL,
    build_analysis_node,
    build_edges,
    build_file_node,
    build_output_file_record,
    get_runtime_context,
    mirror_graph_payload,
    stable_operation_id,
    write_canonical_json,
)


AGE_CODE_MAP = {
    "1": "20-29",
    "2": "30-39",
    "3": "40-49",
    "4": "50-59",
    "5": "60-69",
    "6": "70-79",
}

SEX_CODE_MAP = {
    "1": "M",
    "2": "F",
}


def open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_maybe_gzip(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def normalize_age(raw: str) -> str:
    value = str(raw or "").strip()
    if not value:
        return ""
    if value in AGE_CODE_MAP:
        return AGE_CODE_MAP[value]
    if value in AGE_CODE_MAP.values():
        return value
    compact = value.replace(" ", "")
    if compact in AGE_CODE_MAP.values():
        return compact
    return value


def normalize_sex(raw: str) -> str:
    value = str(raw or "").strip()
    if not value:
        return ""
    return SEX_CODE_MAP.get(value, value)


def derive_subject_id(sample_id: str) -> str:
    value = str(sample_id or "").strip()
    if not value:
        return ""
    parts = value.split("-")
    if len(parts) >= 2:
        return "-".join(parts[:2])
    return ""


def compact_age_comparison_label(age_bin: str, reference_age_bin: str) -> str:
    # Keep the helper name for compatibility with existing callers, but use the
    # readable age-pair identifier directly as comparison_id.
    return expanded_age_comparison_label(age_bin, reference_age_bin)


def expanded_age_comparison_label(age_bin: str, reference_age_bin: str) -> str:
    case_age = str(age_bin or "").strip()
    reference_age = str(reference_age_bin or "").strip()
    if not case_age or not reference_age:
        raise ValueError("age_bin and reference_age_bin must be non-empty")
    return f"{reference_age}_{case_age}"


def parse_gct_header(path: Path) -> tuple[list[str], list[str]]:
    with open_maybe_gzip(path) as handle:
        _version = handle.readline()
        _dims = handle.readline()
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
    if len(header) < 3 or header[0] != "Name":
        raise ValueError("Expected GTEx GCT with Name, Description, and sample columns")
    sample_ids = [str(value).strip() for value in header[2:] if str(value).strip()]
    return header, sample_ids


def build_sample_metadata_rows(
    *,
    sample_rows: list[dict[str, str]],
    subject_rows: list[dict[str, str]],
    gct_sample_ids: list[str],
    tissue_label: str,
    tissue_id: str,
    tissue_column: str | None,
    tissue_value: str | None,
    sample_id_column: str = "SAMPID",
    subject_id_column_sample: str = "SUBJID",
    subject_id_column_subject: str = "SUBJID",
    age_column: str = "AGE",
    sex_column: str = "SEX",
    primary_tissue_column: str = "SMTS",
    detailed_tissue_column: str = "SMTSD",
    age_bins: list[str] | None = None,
) -> list[dict[str, str]]:
    subject_by_id = {
        str(row.get(subject_id_column_subject, "")).strip(): row
        for row in subject_rows
        if str(row.get(subject_id_column_subject, "")).strip()
    }
    gct_sample_id_set = set(gct_sample_ids)
    age_order = list(age_bins or ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79"])

    prepared_meta: list[dict[str, str]] = []
    for row in sample_rows:
        sample_id = str(row.get(sample_id_column, "")).strip()
        if not sample_id or sample_id not in gct_sample_id_set:
            continue
        if tissue_column and tissue_value and str(row.get(tissue_column, "")).strip() != tissue_value:
            continue
        subject_id = str(row.get(subject_id_column_sample, "")).strip()
        if not subject_id:
            subject_id = derive_subject_id(sample_id)
        subject_row = subject_by_id.get(subject_id, {})
        age_bin = normalize_age(str(subject_row.get(age_column, "")))
        sex = normalize_sex(str(subject_row.get(sex_column, "")))
        if age_bin not in age_order:
            continue
        prepared_meta.append(
            {
                "sample_id": sample_id,
                "subject_id": subject_id,
                "age_bin": age_bin,
                "SEX": sex,
                "primary_tissue": str(row.get(primary_tissue_column, "")).strip(),
                "detailed_tissue": str(row.get(detailed_tissue_column, "")).strip(),
                "tissue_id": tissue_id,
                "tissue_label": tissue_label,
            }
        )
    return prepared_meta


def write_filtered_counts(
    *,
    counts_gct: Path,
    sample_columns: list[str],
    sample_index: list[int],
    out_path: Path,
) -> int:
    fieldnames = ["gene_id", "gene_symbol", *sample_columns]
    n_rows = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        with open_maybe_gzip(counts_gct) as in_handle:
            _version = in_handle.readline()
            _dims = in_handle.readline()
            reader = csv.reader(in_handle, delimiter="\t")
            _header = next(reader)
            for row in reader:
                out_row = {
                    "gene_id": str(row[0]).strip(),
                    "gene_symbol": str(row[1]).strip(),
                }
                for out_col, idx in zip(sample_columns, sample_index):
                    out_row[out_col] = row[idx + 2]
                writer.writerow(out_row)
                n_rows += 1
    return n_rows


def build_age_binned_comparisons(
    *,
    prepared_meta: list[dict[str, str]],
    reference_age_bin: str,
    age_bins: list[str],
    min_samples_per_group: int,
) -> tuple[list[dict[str, str]], dict[str, int]]:
    age_counts: dict[str, int] = {}
    for row in prepared_meta:
        age_counts[row["age_bin"]] = age_counts.get(row["age_bin"], 0) + 1

    comparisons: list[dict[str, str]] = []
    for age_bin in age_bins:
        if age_bin == reference_age_bin:
            continue
        if age_counts.get(reference_age_bin, 0) < int(min_samples_per_group):
            continue
        if age_counts.get(age_bin, 0) < int(min_samples_per_group):
            continue
        comparisons.append(
            {
                "comparison_id": compact_age_comparison_label(age_bin, reference_age_bin),
                "gmt_comparison_label": expanded_age_comparison_label(age_bin, reference_age_bin),
                "comparison_kind": "condition_a_vs_b",
                "group_column": "age_bin",
                "group_a": age_bin,
                "group_b": reference_age_bin,
            }
        )
    return comparisons, age_counts


def write_workflow_provenance_graph(
    *,
    workflow_name: str,
    module_name: str,
    output_dir: Path,
    focus_output_path: Path,
    output_paths: list[tuple[Path, str]],
    input_paths: list[tuple[Path, str]],
    parameters: dict[str, Any],
    analysis_description: str | None = None,
    input_overlays: dict[str, dict[str, Any]] | None = None,
) -> Path:
    runtime_ctx = get_runtime_context()
    mirror_local_prefix = runtime_ctx.provenance_mirror_local_prefix if runtime_ctx is not None else None
    mirror_remote_prefix = runtime_ctx.provenance_mirror_remote_prefix if runtime_ctx is not None else None
    # Wrap the caller's {role: overlay} map into the {"inputs": {"role:<role>": ...}}
    # shape that build_file_node / _overlay_for_file expects, so per-input source
    # identifiers (canonical_uri, provider, source) land on the file nodes.
    overlay_arg = {"inputs": {f"role:{role}": ov for role, ov in (input_overlays or {}).items()}}
    input_records = [{"path": str(path), "role": role} for path, role in input_paths]
    input_nodes = [
        build_file_node(
            record,
            overlay_arg,
            mirror_local_prefix=mirror_local_prefix,
            mirror_remote_prefix=mirror_remote_prefix,
        )
        for record in input_records
    ]
    output_records = [
        build_output_file_record(output_dir, {"path": str(path), "role": role})
        for path, role in output_paths
    ]
    output_nodes = [
        build_file_node(
            record,
            {},
            mirror_local_prefix=mirror_local_prefix,
            mirror_remote_prefix=mirror_remote_prefix,
        )
        for record in output_records
    ]
    focus_node = next(node for node in output_nodes if node["name"] == focus_output_path.name)
    extra_output_nodes = [node for node in output_nodes if node["id"] != focus_node["id"]]
    invocation = current_invocation_context()
    command = invocation.get("argv") if invocation else list(sys.argv)
    entrypoint = f"geneset-extractors workflows {workflow_name}"
    operation_id = stable_operation_id(
        workflow_name,
        str(focus_output_path),
        [str(node["id"]) for node in input_nodes],
    )
    operation = build_analysis_node(
        analysis_id=operation_id,
        method=workflow_name,
        name=f"prepare_{focus_output_path.stem}",
        description=analysis_description
        or f"Analysis step that prepares GTEx differential expression results and emits {focus_output_path.name}.",
        parameters=parameters,
        command=command,
        entrypoint=entrypoint,
        repo_url=REPO_URL,
        module=module_name,
        script_url=REPO_URL,
        version=_resolve_git_commit(),
        dcc_url=REPO_URL,
        drc_url=REPO_URL,
    )
    graph_path = output_dir / f"{focus_output_path.stem}.provenance_graph.json"
    payload = mirror_graph_payload(
        {
            focus_output_path.stem: {
                "nodes": input_nodes + [operation] + output_nodes,
                "edges": build_edges(input_nodes, str(operation["id"]), str(focus_node["id"]), extra_output_nodes),
            }
        },
        mirror_local_prefix,
        mirror_remote_prefix,
    )
    write_canonical_json(graph_path, payload)
    return graph_path
