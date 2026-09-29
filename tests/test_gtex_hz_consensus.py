from __future__ import annotations

from geneset_extractors.workflows.gtex_hz_consensus import _groups, _sample_up_sets
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
