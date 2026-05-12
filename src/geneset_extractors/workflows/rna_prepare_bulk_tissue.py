from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import sys
from typing import Iterable

from geneset_extractors.core.metadata import current_invocation_context, input_file_record
from geneset_extractors.core.provenance import (
    REPO_URL,
    build_analysis_node,
    build_edges,
    build_file_node,
    build_output_file_record,
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
    left = str(age_bin or "").strip()
    right = str(reference_age_bin or "").strip()
    if not left or not right:
        raise ValueError("age_bin and reference_age_bin must be non-empty")
    left_decade = left.split("-", 1)[0]
    right_decade = right.split("-", 1)[0]
    return f"age{left_decade}_{right_decade}"


def write_naming_reference(path: Path) -> None:
    text = """# Naming Reference

This prepared bulk RNA bundle uses compact age-bin comparison names.

## Comparison Labels

- `age30_20` means `30-39` vs `20-29`
- `age40_20` means `40-49` vs `20-29`
- `age50_20` means `50-59` vs `20-29`
- `age60_20` means `60-69` vs `20-29`
- `age70_20` means `70-79` vs `20-29`

The suffix `_20` always refers to the reference age bin `20-29`.
"""
    path.write_text(text, encoding="utf-8")


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_gct(path: Path) -> tuple[list[str], list[list[str]]]:
    with open_maybe_gzip(path) as handle:
        _version = handle.readline()
        _dims = handle.readline()
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = [row for row in reader]
    return header, rows


def _append_cli_value(argv: list[str], name: str, value: object) -> None:
    if value is None:
        return
    if isinstance(value, bool):
        argv.extend([f"--{name}", "true" if value else "false"])
        return
    argv.extend([f"--{name}", str(value)])


def _resolve_graph_path(tsv_path: Path) -> Path:
    name = tsv_path.name
    if name.endswith(".tsv.gz"):
        stem = name[:-7]
    elif name.endswith(".tsv"):
        stem = name[:-4]
    else:
        stem = tsv_path.stem
    return tsv_path.with_name(f"{stem}.provenance_graph.json")


def _reconstruct_prepare_command(args) -> list[str]:
    argv = [sys.executable, "-m", "geneset_extractors.cli", "workflows", "rna_prepare_bulk_tissue"]
    _append_cli_value(argv, "counts_gct", args.counts_gct)
    _append_cli_value(argv, "sample_metadata_tsv", args.sample_metadata_tsv)
    _append_cli_value(argv, "subject_metadata_tsv", args.subject_metadata_tsv)
    _append_cli_value(argv, "tissue_label", args.tissue_label)
    _append_cli_value(argv, "out_dir", args.out_dir)
    _append_cli_value(argv, "sample_id_column", args.sample_id_column)
    _append_cli_value(argv, "subject_id_column_sample", args.subject_id_column_sample)
    _append_cli_value(argv, "subject_id_column_subject", args.subject_id_column_subject)
    _append_cli_value(argv, "age_column", args.age_column)
    _append_cli_value(argv, "sex_column", args.sex_column)
    _append_cli_value(argv, "primary_tissue_column", args.primary_tissue_column)
    _append_cli_value(argv, "detailed_tissue_column", args.detailed_tissue_column)
    _append_cli_value(argv, "reference_age_bin", args.reference_age_bin)
    _append_cli_value(argv, "age_bins", args.age_bins)
    _append_cli_value(argv, "min_samples_per_group", args.min_samples_per_group)
    return argv


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_rows = read_tsv(Path(args.sample_metadata_tsv))
    subject_rows = read_tsv(Path(args.subject_metadata_tsv))
    header, gct_rows = parse_gct(Path(args.counts_gct))

    if len(header) < 3 or header[0] != "Name" or header[1] != "Description":
        raise ValueError("Expected GCT with Name, Description, and sample columns")

    subject_by_id = {
        str(row.get(args.subject_id_column_subject, "")).strip(): row
        for row in subject_rows
        if str(row.get(args.subject_id_column_subject, "")).strip()
    }

    gct_sample_ids = [str(value).strip() for value in header[2:] if str(value).strip()]
    gct_sample_id_set = set(gct_sample_ids)
    age_order = [token.strip() for token in str(args.age_bins).split(",") if token.strip()]

    prepared_meta: list[dict[str, str]] = []
    retained_sample_ids: list[str] = []
    for row in sample_rows:
        sample_id = str(row.get(args.sample_id_column, "")).strip()
        if not sample_id or sample_id not in gct_sample_id_set:
            continue
        subject_id = str(row.get(args.subject_id_column_sample, "")).strip()
        if not subject_id:
            subject_id = derive_subject_id(sample_id)
        subject_row = subject_by_id.get(subject_id, {})
        age_bin = normalize_age(str(subject_row.get(args.age_column, "")))
        sex = normalize_sex(str(subject_row.get(args.sex_column, "")))
        if age_bin not in age_order:
            continue
        retained_sample_ids.append(sample_id)
        prepared_meta.append(
            {
                "sample_id": sample_id,
                "subject_id": subject_id,
                "age_bin": age_bin,
                "SEX": sex,
                "primary_tissue": str(row.get(args.primary_tissue_column, "")).strip(),
                "detailed_tissue": str(row.get(args.detailed_tissue_column, "")).strip(),
            }
        )

    retained_sample_id_set = set(retained_sample_ids)
    sample_index = [idx for idx, sample_id in enumerate(gct_sample_ids) if sample_id in retained_sample_id_set]
    sample_columns = [gct_sample_ids[idx] for idx in sample_index]

    counts_rows: list[dict[str, str]] = []
    for row in gct_rows:
        gene_id = str(row[0]).strip()
        gene_symbol = str(row[1]).strip()
        out_row = {"gene_id": gene_id, "gene_symbol": gene_symbol}
        for out_col, idx in zip(sample_columns, sample_index):
            out_row[out_col] = row[idx + 2]
        counts_rows.append(out_row)

    counts_fieldnames = ["gene_id", "gene_symbol", *sample_columns]
    counts_path = out_dir / "tissue_counts.tsv"
    sample_metadata_path = out_dir / "sample_metadata.tsv"
    comparisons_path = out_dir / "comparisons.tsv"
    summary_path = out_dir / "prepare_summary.json"
    naming_reference_path = out_dir / "naming_reference.md"

    write_tsv(counts_path, counts_rows, counts_fieldnames)
    write_tsv(
        sample_metadata_path,
        prepared_meta,
        ["sample_id", "subject_id", "age_bin", "SEX", "primary_tissue", "detailed_tissue"],
    )

    age_counts: dict[str, int] = {}
    for row in prepared_meta:
        age_counts[row["age_bin"]] = age_counts.get(row["age_bin"], 0) + 1

    comparisons: list[dict[str, str]] = []
    reference_age = str(args.reference_age_bin)
    for age_bin in age_order:
        if age_bin == reference_age:
            continue
        if age_counts.get(reference_age, 0) < int(args.min_samples_per_group):
            continue
        if age_counts.get(age_bin, 0) < int(args.min_samples_per_group):
            continue
        comparisons.append(
            {
                "comparison_id": compact_age_comparison_label(age_bin, reference_age),
                "comparison_kind": "condition_a_vs_b",
                "group_column": "age_bin",
                "group_a": age_bin,
                "group_b": reference_age,
            }
        )

    write_tsv(
        comparisons_path,
        comparisons,
        ["comparison_id", "comparison_kind", "group_column", "group_a", "group_b"],
    )

    summary = {
        "workflow": "rna_prepare_bulk_tissue",
        "tissue_label": args.tissue_label,
        "counts_gct": str(Path(args.counts_gct).resolve()),
        "sample_metadata_tsv": str(Path(args.sample_metadata_tsv).resolve()),
        "subject_metadata_tsv": str(Path(args.subject_metadata_tsv).resolve()),
        "n_samples_retained": len(prepared_meta),
        "n_genes_retained": len(counts_rows),
        "reference_age_bin": reference_age,
        "age_bin_counts": {key: age_counts.get(key, 0) for key in age_order},
        "n_comparisons": len(comparisons),
        "comparison_ids": [row["comparison_id"] for row in comparisons],
        "out_dir": str(out_dir),
    }
    write_naming_reference(naming_reference_path)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    input_records = [
        input_file_record(Path(args.counts_gct), "counts_gct"),
        input_file_record(Path(args.sample_metadata_tsv), "sample_metadata_tsv"),
        input_file_record(Path(args.subject_metadata_tsv), "subject_metadata_tsv"),
    ]
    output_records = [
        build_output_file_record(out_dir, {"path": str(counts_path), "role": "counts_tsv"}),
        build_output_file_record(out_dir, {"path": str(sample_metadata_path), "role": "sample_metadata_tsv"}),
        build_output_file_record(out_dir, {"path": str(comparisons_path), "role": "comparisons_tsv"}),
        build_output_file_record(out_dir, {"path": str(summary_path), "role": "prepare_summary"}),
        build_output_file_record(out_dir, {"path": str(naming_reference_path), "role": "naming_reference"}),
    ]
    input_nodes = [build_file_node(record, {}) for record in input_records]
    output_nodes = [build_file_node(record, {}) for record in output_records]
    counts_node = next(node for node in output_nodes if node["name"] == counts_path.name)
    extra_output_nodes = [node for node in output_nodes if node["id"] != counts_node["id"]]
    invocation = current_invocation_context()
    command_argv = [str(token) for token in invocation.get("argv", [])] if invocation else []
    if not command_argv or "rna_prepare_bulk_tissue" not in " ".join(command_argv):
        command_argv = _reconstruct_prepare_command(args)
    operation_id = stable_operation_id(
        "rna_prepare_bulk_tissue",
        str(counts_path),
        [str(node["id"]) for node in input_nodes],
    )
    operation = build_analysis_node(
        analysis_id=operation_id,
        method="rna_prepare_bulk_tissue",
        name=f"prepare_{counts_path.stem}",
        description="Analysis step that prepares a bulk RNA tissue bundle from a raw GCT matrix.",
        parameters={
            "tissue_label": args.tissue_label,
            "sample_id_column": args.sample_id_column,
            "subject_id_column_sample": args.subject_id_column_sample,
            "subject_id_column_subject": args.subject_id_column_subject,
            "age_column": args.age_column,
            "sex_column": args.sex_column,
            "primary_tissue_column": args.primary_tissue_column,
            "detailed_tissue_column": args.detailed_tissue_column,
            "reference_age_bin": args.reference_age_bin,
            "age_bins": [token.strip() for token in str(args.age_bins).split(",") if token.strip()],
            "min_samples_per_group": int(args.min_samples_per_group),
        },
        command=command_argv,
        entrypoint="geneset-extractors workflows rna_prepare_bulk_tissue",
        repo_url=REPO_URL,
        module="geneset_extractors.workflows.rna_prepare_bulk_tissue",
        script_url=REPO_URL,
        version=None,
        dcc_url=REPO_URL,
        drc_url=REPO_URL,
    )
    graph_path = _resolve_graph_path(counts_path)
    write_canonical_json(
        graph_path,
        {
            "tissue_counts": {
                "nodes": input_nodes + [operation] + output_nodes,
                "edges": build_edges(input_nodes, str(operation["id"]), str(counts_node["id"]), extra_output_nodes),
            }
        },
    )
    summary["counts_provenance_graph_path"] = str(graph_path)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary
