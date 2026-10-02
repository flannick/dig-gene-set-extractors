from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import h5py
import numpy as np
import pytest

from geneset_extractors.workflows.lincs_l1000_cp import (
    inspect_gctx,
    merge_partitions,
    plan_cell_time_partitions,
    resolve_last_occurrences,
    run,
)


def _gctx(path: Path) -> None:
    genes = np.asarray([f"GENE{index:03d}".encode() for index in range(500)])
    signatures = np.asarray([b"fixture_a", b"fixture_b"])
    matrix = np.vstack((np.arange(500, 0, -1), np.arange(500)))
    with h5py.File(path, "w") as handle:
        handle.create_dataset("0/DATA/0/matrix", data=matrix)
        handle.create_dataset("0/META/ROW/id", data=genes)
        handle.create_dataset("0/META/COL/lincs_id", data=signatures)
        handle.create_dataset("0/META/COL/cell_line", data=np.asarray([b"CELL_A", b"CELL_B"]))
        handle.create_dataset("0/META/COL/pert_time", data=np.asarray([b"24 h", b"6 h"]))


def _args(gctx_path: Path, out_dir: Path, start: int = 0, end: int | None = None) -> Namespace:
    return Namespace(gctx_path=str(gctx_path), out_dir=str(out_dir), top_n=250, block_size=1, start_index=start, end_index=end)


def test_validates_and_exports_gctx_partitions(tmp_path: Path) -> None:
    gctx_path = tmp_path / "fixture.gctx"
    _gctx(gctx_path)
    genes, signatures, shape = inspect_gctx(gctx_path)
    assert (len(genes), signatures, shape) == (500, ["fixture_a", "fixture_b"], (2, 500))
    result = run(_args(gctx_path, tmp_path / "part_a", 0, 1))
    assert result["n_sets"] == 2
    lines = (tmp_path / "part_a/l1000_cp.gmt").read_text(encoding="utf-8").splitlines()
    assert [line.split("\t", 1)[0] for line in lines] == ["fixture_a up", "fixture_a down"]
    assert all(len(line.split("\t")) == 252 for line in lines)


def test_merge_requires_contiguous_partitions(tmp_path: Path) -> None:
    gctx_path = tmp_path / "fixture.gctx"
    _gctx(gctx_path)
    first, second = tmp_path / "first", tmp_path / "second"
    run(_args(gctx_path, first, 0, 1))
    run(_args(gctx_path, second, 1, 2))
    merged = merge_partitions([second, first], tmp_path / "merged")
    assert merged == {"n_partitions": 2, "n_terms": 4, "end_index": 2}
    assert len((tmp_path / "merged/l1000_cp.gmt").read_text(encoding="utf-8").splitlines()) == 4


def test_duplicate_lincs_ids_retain_the_last_gctx_column(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.gctx"
    genes = np.asarray([f"GENE{index:03d}".encode() for index in range(500)])
    with h5py.File(path, "w") as handle:
        handle.create_dataset("0/DATA/0/matrix", data=np.vstack((np.arange(500, 0, -1), np.arange(500))))
        handle.create_dataset("0/META/ROW/id", data=genes)
        handle.create_dataset("0/META/COL/lincs_id", data=np.asarray([b"duplicate", b"duplicate"]))
        handle.create_dataset("0/META/COL/cell_line", data=np.asarray([b"CELL_A", b"CELL_A"]))
        handle.create_dataset("0/META/COL/pert_time", data=np.asarray([b"24 h", b"24 h"]))
    _, signatures, _ = inspect_gctx(path)
    retained, groups = resolve_last_occurrences(signatures)
    assert retained == [1]
    assert groups == {"duplicate": [0, 1]}
    run(_args(path, tmp_path / "out"))
    lines = (tmp_path / "out/l1000_cp.gmt").read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2
    assert lines[0].split("\t")[2] == "GENE499"
    summary = __import__("json").loads((tmp_path / "out/lincs_l1000_cp_partition.json").read_text(encoding="utf-8"))["duplicate_resolution"]
    assert summary["unique_lincs_id"] == 1
    assert summary["differing_duplicate_vector_groups"] == 1


def test_plans_cell_line_time_worklists_and_exports_one_task(tmp_path: Path) -> None:
    gctx_path = tmp_path / "fixture.gctx"
    _gctx(gctx_path)
    plan_dir = tmp_path / "plan"
    assert plan_cell_time_partitions(gctx_path, plan_dir, max_signatures_per_task=1) == 2
    rows = list(__import__("csv").DictReader((plan_dir / "task_manifest.tsv").open(encoding="utf-8"), delimiter="\t"))
    assert [(row["cell_line"], row["pert_time"], row["n_signatures"]) for row in rows] == [("CELL_A", "24 h", "1"), ("CELL_B", "6 h", "1")]
    args = _args(gctx_path, tmp_path / "task")
    args.raw_indices_tsv = rows[1]["raw_indices_tsv"]
    result = run(args)
    assert result["n_sets"] == 2
    assert (tmp_path / "task/l1000_cp.gmt").read_text(encoding="utf-8").splitlines()[0].startswith("fixture_b up\t")


def test_rejects_missing_required_dataset(tmp_path: Path) -> None:
    path = tmp_path / "bad.gctx"
    with h5py.File(path, "w") as handle:
        handle.create_dataset("0/DATA/0/matrix", data=np.zeros((1, 1)))
    with pytest.raises(ValueError, match="missing required"):
        inspect_gctx(path)
