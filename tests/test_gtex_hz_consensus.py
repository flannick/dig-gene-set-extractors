from __future__ import annotations

from geneset_extractors.workflows.gtex_hz_consensus import _ecdf_up_gene_indices, _groups, _sample_up_sets
import numpy as np


def test_groups_filter_to_expression_samples_and_require_three_members() -> None:
    rows = [
        {"sample_id": f"S{i}", "detailed_tissue": "Brain - Cortex", "SEX": "2", "age_bin": "40-49"}
        for i in range(1, 4)
    ] + [{"sample_id": "S4", "detailed_tissue": "Liver", "SEX": "1", "age_bin": "50-59"}]
    assert _groups(rows, {"S1", "S2", "S3"}, 3) == {("Brain - Cortex", "Female", "40-49"): ["S1", "S2", "S3"]}


def test_sample_up_calls_are_deterministic() -> None:
    calls = _sample_up_sets(["A", "B", "C"], np.array([[1, 2], [4, 2], [8, 2]], dtype=float), 0.95)
    assert calls["0"] == {"C"}
    assert calls["1"] == {"C"}


def test_rank_streaming_up_calls_match_full_signature_transform() -> None:
    symbols = ["A", "B", "C", "D", "E"]
    values = np.array([[1.0], [4.0], [2.0], [8.0], [3.0]])
    full_transform = _sample_up_sets(symbols, values, 0.80)["0"]
    streamed = {symbols[index] for index in _ecdf_up_gene_indices(values[:, 0], 0.80)}
    assert streamed == full_transform
