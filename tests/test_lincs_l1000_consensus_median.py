from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import h5py
import numpy as np

from geneset_extractors.workflows.lincs_l1000_consensus_median import (
    eligible_perturbagens,
    group_signature_indices,
    merge_partitions,
    rank_consensus,
    run,
)


def _gctx(path: Path) -> None:
    genes = np.asarray([f"GENE{index:03d}".encode() for index in range(500)])
    pert_names = np.asarray([b"drug_a"] * 10 + [b"drug_b"] * 9 + [b"drug_c"] * 10)
    rows = []
    for index in range(len(pert_names)):
        values = np.arange(500, dtype=float) + index
        rows.append(values if index < 10 else -values)
    with h5py.File(path, "w") as handle:
        handle.create_dataset("0/DATA/0/matrix", data=np.asarray(rows))
        handle.create_dataset("0/META/ROW/id", data=genes)
        handle.create_dataset("0/META/COL/pert_name", data=pert_names)


def _args(gctx_path: Path, out_dir: Path, index: int = 0, count: int = 1) -> Namespace:
    return Namespace(gctx_path=str(gctx_path), out_dir=str(out_dir), top_n=200, min_signatures=10, partition_index=index, partition_count=count)


def test_threshold_and_deterministic_ranking() -> None:
    groups = group_signature_indices(["a"] * 10 + ["b"] * 9)
    assert eligible_perturbagens(groups, 10) == ["a"]
    genes = ["B", "A", "C", "D"]
    up, down = rank_consensus(genes, np.asarray([1.0, 1.0, 0.0, -1.0]), 2)
    assert up == ["A", "B"]
    assert down == ["C", "D"]


def test_consensus_partitions_and_merge(tmp_path: Path) -> None:
    gctx_path = tmp_path / "fixture.gctx"
    _gctx(gctx_path)
    first, second = tmp_path / "first", tmp_path / "second"
    assert run(_args(gctx_path, first, 0, 2))["n_sets"] == 2
    assert run(_args(gctx_path, second, 1, 2))["n_sets"] == 2
    first_terms = [line.split("\t", 1)[0] for line in (first / "lincs_l1000_consensus_median.gmt").read_text(encoding="utf-8").splitlines()]
    assert first_terms == ["drug_a_up", "drug_a_dn"]
    merged = merge_partitions([second, first], tmp_path / "merged")
    assert merged == {"n_partitions": 2, "n_terms": 4}
