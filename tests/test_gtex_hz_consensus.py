from __future__ import annotations

from geneset_extractors.workflows.gtex_hz_consensus import (
    _gene_ecdf_memmap,
    _groups,
    _modified_row_zscore_memmap,
    _quantile_normalize_by_sample,
    _quantile_normalize_memmap,
    _sample_up_sets,
    _sample_up_support_from_memmap,
)
import numpy as np


def test_groups_filter_to_expression_samples_and_require_three_members() -> None:
    rows = [
        {"sample_id": f"S{i}", "detailed_tissue": "Brain - Cortex", "SEX": "2", "age_bin": "40-49"}
        for i in range(1, 4)
    ] + [{"sample_id": "S4", "detailed_tissue": "Liver", "SEX": "1", "age_bin": "50-59"}]
    assert _groups(rows, {"S1", "S2", "S3"}, 3) == {("Brain - Cortex", "Female", "40-49"): ["S1", "S2", "S3"]}


def test_sample_up_calls_are_deterministic() -> None:
    values = np.array([[1, 2], [4, 2], [8, 2]], dtype=float)
    calls = _sample_up_sets(["A", "B", "C"], values, 0.95)
    assert calls == _sample_up_sets(["A", "B", "C"], values, 0.95)
    assert all(isinstance(genes, set) for genes in calls.values())


def test_reference_quantile_normalization_assigns_ties_to_first_rank() -> None:
    values = np.array([[1.0, 1.0], [1.0, 3.0], [4.0, 5.0]])
    normalized = _quantile_normalize_by_sample(values)
    assert np.array_equal(normalized[:, 0], np.array([1.0, 1.0, 4.5]))
    assert np.array_equal(normalized[:, 1], np.array([1.0, 2.0, 4.5]))


def test_disk_backed_transform_matches_complete_two_stage_reference(tmp_path) -> None:
    symbols = ["A", "B", "C", "D", "E"]
    values = np.array([
        [1.0, 5.0, 2.0, 7.0],
        [4.0, 1.0, 9.0, 3.0],
        [2.0, 6.0, 1.0, 8.0],
        [8.0, 2.0, 7.0, 1.0],
        [3.0, 9.0, 4.0, 2.0],
    ])
    expected = _sample_up_sets(symbols, values, 0.80)
    expression = np.memmap(tmp_path / "expression.mmap", mode="w+", dtype=np.float32, shape=values.shape)
    expression[:, :] = values
    _quantile_normalize_memmap(expression)
    _modified_row_zscore_memmap(expression)
    _, stage_one_cutoff = _gene_ecdf_memmap(expression, values.shape[1], 1, 0.80)
    support, _, _, _ = _sample_up_support_from_memmap(
        expression, stage_one_cutoff, np.arange(values.shape[1]), values.shape[1],
    )
    observed = {
        str(sample_index): {symbols[gene_index] for gene_index in np.flatnonzero(support[sample_index, :])}
        for sample_index in range(values.shape[1])
    }
    assert observed == expected


def test_two_stage_processing_does_not_reduce_to_raw_within_sample_rank() -> None:
    symbols = ["A", "B", "C", "D"]
    values = np.array([
        [100.0, 1.0, 1.0, 1.0],
        [90.0, 2.0, 2.0, 2.0],
        [80.0, 3.0, 3.0, 3.0],
        [70.0, 4.0, 4.0, 4.0],
    ])
    calls = _sample_up_sets(symbols, values, 0.95)
    raw_top = {"A"}
    assert calls["0"] != raw_top
