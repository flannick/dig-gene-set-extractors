"""GTEx V8 tissue, sex, and age consensus gene-set reconstruction.

This workflow owns the biological processing described in the transferred
legacy handoff.  The GTEx wrapper only selects HZ2 and invokes this module.
"""
from __future__ import annotations

import csv
import gzip
import tempfile
import re
from collections import OrderedDict, defaultdict
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
    """Match the reference's column quantile normalization, including ties."""
    sorted_values = np.sort(values, axis=0)
    mean_by_rank = sorted_values.mean(axis=1)
    normalized = np.empty_like(values, dtype=float)
    for column in range(values.shape[1]):
        # The reference assigns every tied value to its *first* sorted rank.
        first_rank = np.searchsorted(sorted_values[:, column], values[:, column], side="left")
        normalized[:, column] = mean_by_rank[first_rank]
    return normalized


def _sample_up_sets(gene_symbols: list[str], values: np.ndarray, cutoff: float) -> dict[str, set[str]]:
    """Return Up calls from the complete historical-style Harmonizome transform.

    ``values`` is gene by sample.  This small in-memory implementation is
    intentionally kept as the reference for the disk-backed production path.
    """
    standardized = _harmonizome_standardize(values)
    up: dict[str, set[str]] = {}
    for column in range(standardized.shape[1]):
        up[str(column)] = {gene_symbols[index] for index in np.flatnonzero(standardized[:, column] >= cutoff)}
    return up


def _ecdf(values: np.ndarray) -> np.ndarray:
    """Right-continuous empirical CDF values, preserving equal-value ties."""
    ordered = np.sort(values, kind="mergesort")
    return np.searchsorted(ordered, values, side="right").astype(np.float64) / len(values)


def _harmonizome_standardize(values: np.ndarray) -> np.ndarray:
    """Reference implementation of the patched Harmonizome transformations."""
    imputed = np.asarray(values, dtype=float).copy()
    imputed[imputed == 0] = np.nan
    row_means = np.nanmean(imputed, axis=1)
    row_means[~np.isfinite(row_means)] = 0.0
    missing = np.where(~np.isfinite(imputed))
    imputed[missing] = row_means[missing[0]]
    normalized = _quantile_normalize_by_sample(np.log2(imputed + 1.0))
    median = np.median(normalized, axis=1, keepdims=True)
    mad = np.median(np.abs(normalized - median), axis=1, keepdims=True)
    meanad = np.mean(np.abs(normalized - median), axis=1, keepdims=True)
    robust_z = np.zeros_like(normalized, dtype=float)
    use_mad = mad != 0
    robust_z[use_mad[:, 0], :] = 0.6745 * (normalized[use_mad[:, 0], :] - median[use_mad[:, 0], :]) / mad[use_mad[:, 0], :]
    use_fallback = ~use_mad[:, 0]
    denominator = 1.253314 * meanad[use_fallback, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        robust_z[use_fallback, :] = (normalized[use_fallback, :] - median[use_fallback, :]) / denominator
    robust_z[~np.isfinite(robust_z)] = 0.0
    gene_ecdf = np.empty_like(robust_z, dtype=np.float64)
    for gene_index in range(robust_z.shape[0]):
        ecdf = _ecdf(robust_z[gene_index, :])
        gene_ecdf[gene_index, :] = 2.0 * (ecdf - ecdf.mean())
    global_ecdf = _ecdf(gene_ecdf.reshape(-1))
    return (2.0 * (global_ecdf - global_ecdf.mean())).reshape(gene_ecdf.shape)


def _gct_row_count(expression_gct: Path) -> int:
    with _open_text(expression_gct) as handle:
        handle.readline()
        dimensions = handle.readline().strip().split("\t")
    try:
        return int(dimensions[0])
    except (IndexError, ValueError) as exc:
        raise ValueError("Expected a GCT dimensions line with a row count") from exc


def _write_gene_major_expression(
    expression_gct: Path,
    temp_dir: Path,
) -> tuple[list[str], list[str], np.memmap]:
    """Materialize a disk-backed, gene-major float32 expression matrix.

    The V8 TPM GCT is too large for nested Python lists or multiple dense
    in-memory arrays.  Gene-major layout permits efficient GCT ingestion; a
    The gene-major layout supports the required gene-wise robust scaling and
    first ECDF stage without putting the full matrix in memory.
    """
    _, all_ids = parse_gct_header(expression_gct)
    max_rows = _gct_row_count(expression_gct)
    nonzero_by_sample = np.zeros(len(all_ids), dtype=np.int64)
    with _open_text(expression_gct) as handle:
        handle.readline(); handle.readline()
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            for index in range(len(all_ids)):
                try:
                    nonzero_by_sample[index] += float(row[index + 2]) != 0.0
                except (IndexError, ValueError):
                    pass
    keep_samples = nonzero_by_sample >= int(0.05 * max_rows)
    sample_ids = [sample_id for sample_id, keep in zip(all_ids, keep_samples) if keep]
    if not sample_ids:
        raise ValueError("The Harmonizome 5% non-missing filter removed every expression sample")
    sample_positions = [index + 2 for index, keep in enumerate(keep_samples) if keep]
    gene_major = np.memmap(temp_dir / "source_expression.float32.mmap", mode="w+", dtype=np.float32, shape=(max_rows, len(sample_ids)))
    source_symbols: list[str] = []
    retained_rows = 0
    min_row_nonzero = int(0.05 * len(sample_ids))
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
            for index, position in enumerate(sample_positions):
                try:
                    values[index] = max(float(row[position] or 0.0), 0.0)
                except (IndexError, ValueError):
                    pass
            if np.count_nonzero(values) < min_row_nonzero:
                continue
            missing = values == 0
            if np.any(missing):
                row_mean = float(values[~missing].mean()) if np.any(~missing) else 0.0
                values[missing] = row_mean
            gene_major[retained_rows, :] = values
            source_symbols.append(symbol)
            retained_rows += 1
    if not source_symbols:
        raise ValueError("No usable gene-symbol rows were found in the expression GCT")
    expression = np.memmap(temp_dir / "expression.float32.mmap", mode="w+", dtype=np.float32, shape=(retained_rows, len(sample_ids)))
    expression[:, :] = gene_major[:retained_rows, :]
    expression.flush()
    del gene_major
    return sample_ids, source_symbols, expression


def _quantile_normalize_memmap(expression: np.memmap) -> None:
    """Quantile-normalize sample columns with only one column in memory."""
    n_genes, n_samples = expression.shape
    mean_by_rank = np.zeros(n_genes, dtype=np.float64)
    for sample_index in range(n_samples):
        values = np.log2(np.maximum(np.asarray(expression[:, sample_index], dtype=np.float64), 0.0) + 1.0)
        expression[:, sample_index] = values
        order = np.argsort(values, kind="mergesort")
        mean_by_rank += values[order]
    mean_by_rank /= n_samples
    for sample_index in range(n_samples):
        values = np.asarray(expression[:, sample_index], dtype=np.float64)
        first_rank = np.searchsorted(np.sort(values, kind="mergesort"), values, side="left")
        expression[:, sample_index] = mean_by_rank[first_rank].astype(np.float32)
    expression.flush()


def _modified_row_zscore_memmap(expression: np.memmap) -> None:
    """Apply the reference modified row z-score in place."""
    n_genes, _ = expression.shape
    for gene_index in range(n_genes):
        values = np.asarray(expression[gene_index, :], dtype=np.float64)
        median = np.median(values)
        deviations = np.abs(values - median)
        mad = np.median(deviations)
        if mad != 0:
            zscore = 0.6745 * (values - median) / mad
        else:
            denominator = 1.253314 * float(np.mean(deviations))
            with np.errstate(divide="ignore", invalid="ignore"):
                zscore = (values - median) / denominator
        zscore[~np.isfinite(zscore)] = 0.0
        expression[gene_index, :] = zscore
    expression.flush()


def _merge_duplicate_symbols_memmap(
    expression: np.memmap, source_symbols: list[str], temp_dir: Path,
) -> tuple[list[str], np.memmap]:
    """Match reference mapping/merge behavior when source symbols are labels."""
    rows_by_symbol: OrderedDict[str, list[int]] = OrderedDict()
    for row_index, symbol in enumerate(source_symbols):
        rows_by_symbol.setdefault(symbol, []).append(row_index)
    symbols = list(rows_by_symbol)
    merged = np.memmap(temp_dir / "mapped_merged.float32.mmap", mode="w+", dtype=np.float32, shape=(len(symbols), expression.shape[1]))
    for merged_index, source_rows in enumerate(rows_by_symbol.values()):
        if len(source_rows) == 1:
            merged[merged_index, :] = expression[source_rows[0], :]
        else:
            merged[merged_index, :] = np.mean(expression[source_rows, :], axis=0, dtype=np.float64)
    merged.flush()
    return symbols, merged


def _gene_ecdf_memmap(
    expression: np.memmap, sample_per_gene: int, seed: int, final_up_cutoff: float,
) -> tuple[float, float]:
    """Apply the reference first ECDF stage in place.

    The reference uses a deterministic, per-gene random sample to estimate the
    final global-ECDF cutoff instead of sorting the entire V8-scale matrix.
    """
    n_genes, n_samples = expression.shape
    sampled_count = min(sample_per_gene, n_samples)
    sampled = np.empty((n_genes, sampled_count), dtype=np.float32)
    rng = np.random.default_rng(seed)
    for gene_index in range(n_genes):
        values = np.asarray(expression[gene_index, :], dtype=np.float64)
        ordered = np.sort(values, kind="mergesort")
        ranks = np.searchsorted(ordered, values, side="right")
        ecdf = ranks.astype(np.float32) / n_samples
        stage_one = 2.0 * (ecdf - ecdf.mean(dtype=np.float64))
        expression[gene_index, :] = stage_one
        sampled[gene_index, :] = stage_one[rng.choice(n_samples, size=sampled_count, replace=False)]
    expression.flush()
    global_sample = sampled.reshape(-1)
    global_sample.sort()
    global_ranks = np.searchsorted(global_sample, global_sample, side="right")
    global_ecdf_mean = float(np.mean(global_ranks / len(global_sample)))
    target_ecdf = min(1.0, global_ecdf_mean + 0.5 * final_up_cutoff)
    quantile_index = min(len(global_sample) - 1, max(0, int(np.ceil(target_ecdf * len(global_sample))) - 1))
    return global_ecdf_mean, float(global_sample[quantile_index])


def _sample_up_support_from_memmap(
    expression: np.memmap,
    stage_one_cutoff: float,
    group_indices: np.ndarray,
    n_groups: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Aggregate reference Up fractions and mean stage-one scores by group."""
    n_genes, n_samples = expression.shape
    support_counts = np.zeros((n_groups, n_genes), dtype=np.uint16)
    score_sums = np.zeros((n_groups, n_genes), dtype=np.float32)
    sample_up_counts = np.zeros(n_samples, dtype=np.int32)
    unique_up_genes = 0
    for gene_index in range(n_genes):
        scores = np.asarray(expression[gene_index, :], dtype=np.float64)
        is_up = scores >= stage_one_cutoff
        sample_up_counts += is_up
        if np.any(is_up):
            unique_up_genes += 1
            valid_up = is_up & (group_indices >= 0)
            support_counts[:, gene_index] = np.bincount(group_indices[valid_up], minlength=n_groups)
        valid_scores = group_indices >= 0
        score_sums[:, gene_index] = np.bincount(
            group_indices[valid_scores], weights=scores[valid_scores], minlength=n_groups,
        )
    return support_counts, score_sums, sample_up_counts, unique_up_genes


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
    with tempfile.TemporaryDirectory(prefix="gtex_hz_consensus_", dir=out_dir) as temp_name:
        sample_ids, source_symbols, expression = _write_gene_major_expression(expression_gct, Path(temp_name))
        _quantile_normalize_memmap(expression)
        _modified_row_zscore_memmap(expression)
        symbols, merged_expression = _merge_duplicate_symbols_memmap(expression, source_symbols, Path(temp_name))
        del expression
        global_ecdf_mean, stage_one_cutoff = _gene_ecdf_memmap(
            merged_expression, int(getattr(args, "global_quantile_sample_per_gene", 500)),
            int(getattr(args, "random_seed", 1)), float(args.up_cutoff),
        )
        group_indices = np.asarray([group_for_sample.get(sample_id, -1) for sample_id in sample_ids], dtype=np.intp)
        group_sample_counts = np.bincount(group_indices[group_indices >= 0], minlength=len(ordered_groups))
        if np.any(group_sample_counts == 0):
            raise ValueError("The Harmonizome 5% sample filter removed every sample from one or more tissue-sex-age groups")
        support_counts, score_sums, sample_up_counts, unique_up_genes = _sample_up_support_from_memmap(
            merged_expression, stage_one_cutoff, group_indices, len(ordered_groups),
        )
        del merged_expression
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
                n_members = int(group_sample_counts[group_index])
                eligible = np.flatnonzero(counts.astype(np.float64) / n_members >= float(args.support_fraction))
                mean_scores = score_sums[group_index, :] / n_members
                ranked = sorted(eligible, key=lambda index: (-(int(counts[index]) / n_members), -float(mean_scores[index]), symbols[index]))
                term = f"{tissue} {sex} {age} Up"
                geneset_name = _consensus_geneset_name(tissue, sex, age)
                genes = [symbols[index] for index in ranked[:int(args.top_n)]]
                gmt_sets.append((geneset_name, genes))
                for index in ranked:
                    row = [geneset_name, symbols[index], int(counts[index]), int(counts[index]) / n_members]
                    support_writer.writerow([term, *row[1:]])
                    full_writer.writerow(row)
                for index in ranked[:int(args.top_n)]:
                    selected_writer.writerow([geneset_name, symbols[index], int(counts[index]), int(counts[index]) / n_members])
        write_gmt(gmt_sets, gmt_path)
    graph = write_workflow_provenance_graph(
        workflow_name="gtex_hz_consensus", module_name=__name__, output_dir=out_dir,
        focus_output_path=gmt_path, output_paths=[(gmt_path, "gmt"), (support_path, "gene_support")],
        input_paths=[(expression_gct, "gtex_v8_tpm"), (sample_attributes, "sample_attributes_tsv_v8"), (subject_phenotypes, "subject_phenotypes_tsv_v8")],
        parameters={"grouping": "SMTSD x SEX x AGE; n_samples >= %d" % int(args.min_samples_per_group), "sample_signature": "5%% non-missing filtering; zero-to-row-mean imputation; log2(TPM+1); sample-column quantile normalization; modified row z-score; duplicate-symbol row-mean merge; gene-wise ECDF and deterministic sampled global-ECDF Up >= %.2f" % float(args.up_cutoff), "global_ecdf_sample_per_gene": int(getattr(args, "global_quantile_sample_per_gene", 500)), "global_ecdf_seed": int(getattr(args, "random_seed", 1)), "consensus": "support_fraction >= %.2f; top_n=%d" % (float(args.support_fraction), int(args.top_n)), "historical_drc_aggregation": "reference-aligned Harmonizome reconstruction"},
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
            "parameters": {"support_fraction": float(args.support_fraction), "top_n": int(args.top_n), "up_cutoff": float(args.up_cutoff), "global_quantile_sample_per_gene": int(getattr(args, "global_quantile_sample_per_gene", 500)), "random_seed": int(getattr(args, "random_seed", 1))},
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
        "sample_up_genes": {
            "median": float(np.median(sample_up_counts)),
            "min": int(sample_up_counts.min()),
            "max": int(sample_up_counts.max()),
        },
        "n_unique_genes_ever_up": int(unique_up_genes),
        "global_ecdf_mean_estimate": global_ecdf_mean,
        "stage1_cutoff_for_final_up": stage_one_cutoff,
        "global_quantile_sample_per_gene": int(getattr(args, "global_quantile_sample_per_gene", 500)),
        "random_seed": int(getattr(args, "random_seed", 1)),
    })
    return {"out_dir": str(out_dir), "n_groups": len(groups), "n_samples": len(sample_ids)}
