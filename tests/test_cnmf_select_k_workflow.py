from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from omics2geneset.cli import main
from omics2geneset.workflows.cnmf_select_k import run as run_cnmf_select_k


def _write_stats_npz(base_dir: Path, name: str) -> Path:
    run_dir = base_dir / name
    run_dir.mkdir(parents=True, exist_ok=True)
    stats_path = run_dir / f"{name}.k_selection_stats.df.npz"
    data = np.array(
        [
            [10.0, 0.92, 1.00],
            [15.0, 0.98, 1.10],
            [20.0, 0.96, 1.25],
            [25.0, 0.93, 1.45],
        ],
        dtype=float,
    )
    columns = np.array(["k", "stability", "prediction_error"], dtype=object)
    index = np.array([0, 1, 2, 3], dtype=int)
    np.savez(stats_path, data=data, columns=columns, index=index)
    return stats_path


class _Args:
    def __init__(self, cnmf_output_dir: Path, name: str):
        self.cnmf_output_dir = str(cnmf_output_dir)
        self.name = name
        self.stats_tsv = None
        self.strategy = "largest_stable"
        self.stability_frac_of_max = 0.95
        self.min_stability_abs = 0.0
        self.require_local_max = False
        self.fixed_k = None


def test_cnmf_select_k_largest_stable_writes_artifacts(tmp_path: Path, capsys):
    out = tmp_path / "cnmf_out"
    name = "toy"
    _write_stats_npz(out, name)
    args = _Args(out, name)
    result = run_cnmf_select_k(args)
    assert int(result["selected_k"]) == 20

    captured = capsys.readouterr()
    assert captured.out == "20\n"

    run_dir = out / name
    assert (run_dir / "omics2geneset_k_selection.tsv").exists()
    assert (run_dir / "omics2geneset_selected_k.json").exists()

    payload = json.loads((run_dir / "omics2geneset_selected_k.json").read_text(encoding="utf-8"))
    assert payload["selected_k"] == 20
    assert payload["strategy"] == "largest_stable"
    assert "selection_summary" in payload


def test_cnmf_select_k_respects_stability_frac_of_max(tmp_path: Path, capsys):
    out = tmp_path / "cnmf_out"
    name = "toy"
    _write_stats_npz(out, name)
    args = _Args(out, name)
    args.stability_frac_of_max = 0.99
    result = run_cnmf_select_k(args)
    assert int(result["selected_k"]) == 15

    captured = capsys.readouterr()
    assert captured.out == "15\n"


def test_cnmf_select_k_manual_fixed_k(tmp_path: Path, capsys):
    out = tmp_path / "cnmf_out"
    name = "toy"
    _write_stats_npz(out, name)
    args = _Args(out, name)
    args.strategy = "manual"
    args.fixed_k = 42
    result = run_cnmf_select_k(args)
    assert int(result["selected_k"]) == 42

    captured = capsys.readouterr()
    assert captured.out == "42\n"

    payload = json.loads((out / name / "omics2geneset_selected_k.json").read_text(encoding="utf-8"))
    assert payload["strategy"] == "manual"
    assert payload["thresholds"]["fixed_k"] == 42


def test_cnmf_select_k_cli_prints_integer_only(tmp_path: Path, capsys):
    out = tmp_path / "cnmf_out"
    name = "toy"
    _write_stats_npz(out, name)
    code = main(
        [
            "workflows",
            "cnmf_select_k",
            "--cnmf_output_dir",
            str(out),
            "--name",
            name,
        ]
    )
    assert code == 0
    captured = capsys.readouterr()
    assert captured.out == "20\n"
