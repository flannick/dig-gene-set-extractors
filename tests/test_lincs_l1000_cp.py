from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import h5py
import numpy as np
import pytest

from geneset_extractors.workflows.lincs_l1000_cp import inspect_gctx, merge_partitions, resolve_last_occurrences, run


def _gctx(path: Path) -> None:
    genes = np.asarray([f"GENE{index:03d}".encode() for index in range(500)])
    signatures = np.asarray([b"fixture_a", b"fixture_b"])
    matrix = np.vstack((np.arange(500, 0, -1), np.arange(500)))
    with h5py.File(path, "w") as handle:
        handle.create_dataset("0/DATA/0/matrix", data=matrix)
        handle.create_dataset("0/META/ROW/id", data=genes)
        handle.create_dataset("0/META/COL/lincs_id", data=signatures)


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


def test_rejects_missing_required_dataset(tmp_path: Path) -> None:
    path = tmp_path / "bad.gctx"
    with h5py.File(path, "w") as handle:
        handle.create_dataset("0/DATA/0/matrix", data=np.zeros((1, 1)))
    with pytest.raises(ValueError, match="missing required"):
        inspect_gctx(path)
