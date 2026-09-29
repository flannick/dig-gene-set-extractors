"""GTEx V8 tissue, sex, and age consensus gene-set reconstruction.

This workflow owns the biological processing described in the transferred
legacy handoff.  The GTEx wrapper only selects HZ2 and invokes this module.
"""
from __future__ import annotations

import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from geneset_extractors.core.provenance import activate_runtime_context
from geneset_extractors.workflows.gtex_runtime_common import (
    build_sample_metadata_rows,
    parse_gct_header,
    read_tsv,
    write_workflow_provenance_graph,
)


AGE_BINS = ("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")


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


def _read_expression(expression_gct: Path, selected_ids: set[str]) -> tuple[list[str], list[str], np.ndarray]:
    _, all_ids = parse_gct_header(expression_gct)
    sample_ids = [sample_id for sample_id in all_ids if sample_id in selected_ids]
    if not sample_ids:
        raise ValueError("No metadata-selected samples occur in the expression GCT")
    sample_positions = {sample_id: index + 2 for index, sample_id in enumerate(all_ids)}
    rows: dict[str, list[float]] = {}
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
            values = []
            for sample_id in sample_ids:
                try:
                    values.append(float(row[sample_positions[sample_id]] or 0.0))
                except (IndexError, ValueError):
                    values.append(0.0)
            # Gene-symbol duplicate merge retains the higher-mean expression row.
            previous = rows.get(symbol)
            if previous is None or sum(values) > sum(previous):
                rows[symbol] = values
    if not rows:
        raise ValueError("No usable gene-symbol rows were found in the expression GCT")
    symbols = sorted(rows)
    return sample_ids, symbols, np.asarray([rows[symbol] for symbol in symbols], dtype=float)


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
    selected_ids = {sample_id for members in groups.values() for sample_id in members}
    sample_ids, symbols, expression = _read_expression(expression_gct, selected_ids)
    calls_by_column = _sample_up_sets(symbols, expression, float(args.up_cutoff))
    calls = {sample_id: calls_by_column[str(index)] for index, sample_id in enumerate(sample_ids)}

    gmt_path = out_dir / "genesets.gmt"
    support_path = out_dir / "gene_support.tsv"
    support_rows: list[tuple[str, str, int, float]] = []
    with gmt_path.open("w", encoding="utf-8", newline="\n") as handle:
        for (tissue, sex, age), members in sorted(groups.items(), key=lambda item: f"{item[0][0]} {item[0][1]} {item[0][2]} Up"):
            counts: Counter[str] = Counter(gene for sample_id in members for gene in calls[sample_id])
            ranked = [(gene, count, count / len(members)) for gene, count in counts.items() if count / len(members) >= float(args.support_fraction)]
            ranked.sort(key=lambda item: (-item[2], -item[1], item[0]))
            term = f"{tissue} {sex} {age} Up"
            genes = [gene for gene, _, _ in ranked[:int(args.top_n)]]
            handle.write("\t".join([term, str(args.description), *genes]) + "\n")
            support_rows.extend((term, gene, count, fraction) for gene, count, fraction in ranked)
    with support_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["term", "gene_symbol", "support_count", "support_fraction"])
        writer.writerows(support_rows)
    graph = write_workflow_provenance_graph(
        workflow_name="gtex_hz_consensus", module_name=__name__, output_dir=out_dir,
        focus_output_path=gmt_path, output_paths=[(gmt_path, "gmt"), (support_path, "gene_support")],
        input_paths=[(expression_gct, "gtex_v8_tpm"), (sample_attributes, "sample_attributes_tsv_v8"), (subject_phenotypes, "subject_phenotypes_tsv_v8")],
        parameters={"grouping": "SMTSD x SEX x AGE; n_samples >= %d" % int(args.min_samples_per_group), "sample_signature": "log2(TPM+1), sample quantile normalization, robust median/MAD z-score with mean-absolute-deviation fallback, ECDF Up >= %.2f" % float(args.up_cutoff), "consensus": "support_fraction >= %.2f; top_n=%d" % (float(args.support_fraction), int(args.top_n)), "historical_drc_aggregation": "unavailable; scientifically comparable reconstruction, not set-equivalent"},
    )
    graph.replace(out_dir / "geneset.provenance.legacy.json")
    return {"out_dir": str(out_dir), "n_groups": len(groups), "n_samples": len(sample_ids)}
