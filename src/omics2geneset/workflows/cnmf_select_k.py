from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np


REQUIRED_COLUMNS = ("k", "stability", "prediction_error")


def _as_text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _float_or_none(value: Any) -> float | None:
    if value is None:
        return None
    text = _as_text(value).strip()
    if not text:
        return None
    try:
        out = float(text)
    except ValueError:
        return None
    if out != out:
        return None
    return float(out)


def _load_rows_from_stats_tsv(path: Path) -> list[dict[str, float]]:
    with path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"stats_tsv has no header: {path}")
        field_set = set(str(x) for x in reader.fieldnames)
        missing = [col for col in REQUIRED_COLUMNS if col not in field_set]
        if missing:
            raise ValueError(
                "stats_tsv missing required columns: "
                + ", ".join(missing)
                + ". Required columns are: "
                + ", ".join(REQUIRED_COLUMNS)
            )
        rows: list[dict[str, float]] = []
        for row in reader:
            parsed: dict[str, float] = {}
            skip_row = False
            for key, raw in row.items():
                if key is None:
                    continue
                val = _float_or_none(raw)
                if key in REQUIRED_COLUMNS and val is None:
                    skip_row = True
                    break
                if val is not None:
                    parsed[str(key)] = float(val)
            if skip_row:
                continue
            rows.append(parsed)
    return rows


def _rows_from_data_columns(data: np.ndarray, columns: list[str]) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    if data.ndim != 2:
        raise ValueError("df.npz `data` array must be 2D when `columns` are provided.")
    if data.shape[1] != len(columns):
        raise ValueError(
            f"df.npz mismatch: data has {data.shape[1]} columns but metadata has {len(columns)} column names."
        )
    for row in data:
        parsed: dict[str, float] = {}
        for key, raw in zip(columns, row):
            val = _float_or_none(raw)
            if val is not None:
                parsed[str(key)] = float(val)
        rows.append(parsed)
    return rows


def _rows_from_structured_data(data: np.ndarray) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    if not data.dtype.names:
        raise ValueError("Structured data array expected dtype names but none were present.")
    names = [str(x) for x in data.dtype.names]
    for row in data:
        parsed: dict[str, float] = {}
        for key in names:
            val = _float_or_none(row[key])
            if val is not None:
                parsed[key] = float(val)
        rows.append(parsed)
    return rows


def _rows_from_direct_arrays(npz: np.lib.npyio.NpzFile) -> list[dict[str, float]]:
    keys = list(npz.files)
    if not set(REQUIRED_COLUMNS).issubset(set(keys)):
        raise ValueError("npz did not include direct required arrays for k/stability/prediction_error.")
    lengths = []
    arrays: dict[str, np.ndarray] = {}
    for key in keys:
        arr = np.asarray(npz[key])
        if arr.ndim == 0:
            arr = np.asarray([arr.item()])
        arrays[key] = arr
        lengths.append(arr.shape[0])
    n = min(lengths) if lengths else 0
    rows: list[dict[str, float]] = []
    for i in range(n):
        parsed: dict[str, float] = {}
        for key, arr in arrays.items():
            val = _float_or_none(arr[i])
            if val is not None:
                parsed[str(key)] = float(val)
        rows.append(parsed)
    return rows


def _load_rows_from_df_npz(path: Path) -> list[dict[str, float]]:
    try:
        npz = np.load(path, allow_pickle=True)
    except Exception as exc:
        raise ValueError(f"Failed to read k-selection npz file: {path} ({exc})") from exc

    with npz:
        files = list(npz.files)
        rows: list[dict[str, float]]
        if {"data", "columns"}.issubset(set(files)):
            data = np.asarray(npz["data"])
            cols_raw = np.asarray(npz["columns"])
            columns = [_as_text(x) for x in cols_raw.tolist()]
            rows = _rows_from_data_columns(data, columns)
        elif "data" in files and np.asarray(npz["data"]).dtype.names:
            rows = _rows_from_structured_data(np.asarray(npz["data"]))
        else:
            rows = _rows_from_direct_arrays(npz)

    filtered: list[dict[str, float]] = []
    for row in rows:
        k = _float_or_none(row.get("k"))
        s = _float_or_none(row.get("stability"))
        e = _float_or_none(row.get("prediction_error"))
        if k is None or s is None or e is None:
            continue
        clean = dict(row)
        clean["k"] = float(k)
        clean["stability"] = float(s)
        clean["prediction_error"] = float(e)
        filtered.append(clean)

    if not filtered:
        raise ValueError(
            "Could not parse valid k/stability/prediction_error rows from df.npz. "
            "Re-run `cnmf k_selection_plot` or pass --stats_tsv exported from k-selection stats."
        )
    return filtered


def _is_local_max(sorted_rows: list[dict[str, float]], idx: int) -> bool:
    curr = float(sorted_rows[idx]["stability"])
    if idx > 0 and curr < float(sorted_rows[idx - 1]["stability"]):
        return False
    if idx < len(sorted_rows) - 1 and curr < float(sorted_rows[idx + 1]["stability"]):
        return False
    return True


def _select_max_stability(rows: list[dict[str, float]]) -> tuple[int, str]:
    sorted_rows = sorted(rows, key=lambda r: (float(r["stability"]), float(r["k"])))
    chosen = sorted_rows[-1]
    return int(round(float(chosen["k"]))), "selected max stability (tie-breaker: larger k)"


def _select_largest_stable(
    rows: list[dict[str, float]],
    *,
    stability_frac_of_max: float,
    min_stability_abs: float,
    require_local_max: bool,
) -> tuple[int, str]:
    by_k = sorted(rows, key=lambda r: float(r["k"]))
    s_max = max(float(r["stability"]) for r in by_k)
    threshold = max(float(min_stability_abs), float(s_max) * float(stability_frac_of_max))
    candidates = [r for r in by_k if float(r["stability"]) >= threshold]
    if require_local_max:
        index_by_k = {float(r["k"]): i for i, r in enumerate(by_k)}
        candidates = [r for r in candidates if _is_local_max(by_k, index_by_k[float(r["k"])])]
    if candidates:
        chosen = max(candidates, key=lambda r: float(r["k"]))
        reason = (
            "selected largest k meeting stability thresholds: "
            f"stability >= max({min_stability_abs:.4g}, {stability_frac_of_max:.4g} * s_max={s_max:.4g})"
        )
        if require_local_max:
            reason += " and local-max criterion"
        return int(round(float(chosen["k"]))), reason
    k_fallback, _ = _select_max_stability(by_k)
    return k_fallback, "no k passed largest_stable thresholds; fell back to max stability"


def _write_rows_tsv(path: Path, rows: list[dict[str, float]]) -> list[str]:
    keys = set()
    for row in rows:
        keys.update(str(k) for k in row.keys())
    first_cols = [c for c in REQUIRED_COLUMNS if c in keys]
    rest = sorted(k for k in keys if k not in set(first_cols))
    columns = first_cols + rest
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        for row in sorted(rows, key=lambda r: float(r["k"])):
            out: dict[str, object] = {}
            for col in columns:
                if col in row:
                    out[col] = row[col]
            writer.writerow(out)
    return columns


def run(args) -> dict[str, object]:
    cnmf_output_dir = Path(args.cnmf_output_dir)
    run_dir = cnmf_output_dir / str(args.name)
    run_dir.mkdir(parents=True, exist_ok=True)

    stats_path = run_dir / f"{args.name}.k_selection_stats.df.npz"
    stats_tsv = Path(args.stats_tsv) if getattr(args, "stats_tsv", None) else None

    if stats_tsv is not None:
        rows = _load_rows_from_stats_tsv(stats_tsv)
        source_path = stats_tsv
    else:
        if not stats_path.exists():
            raise ValueError(
                "Could not find cNMF stats file: "
                f"{stats_path}. Re-run `cnmf k_selection_plot` or pass --stats_tsv."
            )
        rows = _load_rows_from_df_npz(stats_path)
        source_path = stats_path

    if not rows:
        raise ValueError("No k-selection rows were available after parsing.")

    strategy = str(args.strategy)
    fixed_k = getattr(args, "fixed_k", None)
    if strategy == "manual":
        if fixed_k is None:
            raise ValueError("--fixed_k is required when --strategy manual")
        selected_k = int(fixed_k)
        reason = f"manual strategy selected fixed_k={selected_k}"
    elif strategy == "max_stability":
        selected_k, reason = _select_max_stability(rows)
    elif strategy == "largest_stable":
        selected_k, reason = _select_largest_stable(
            rows,
            stability_frac_of_max=float(args.stability_frac_of_max),
            min_stability_abs=float(args.min_stability_abs),
            require_local_max=bool(args.require_local_max),
        )
    else:
        raise ValueError(f"Unsupported strategy: {strategy}")

    rows_path = run_dir / "omics2geneset_k_selection.tsv"
    columns = _write_rows_tsv(rows_path, rows)

    payload = {
        "selected_k": int(selected_k),
        "strategy": strategy,
        "thresholds": {
            "stability_frac_of_max": float(args.stability_frac_of_max),
            "min_stability_abs": float(args.min_stability_abs),
            "require_local_max": bool(args.require_local_max),
            "fixed_k": (int(fixed_k) if fixed_k is not None else None),
        },
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "stats_source_path": str(source_path),
        "stats_rows_written": int(len(rows)),
        "stats_columns": columns,
        "selection_summary": reason,
    }
    json_path = run_dir / "omics2geneset_selected_k.json"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    # Must print a single integer only for shell scripting.
    sys.stdout.write(f"{int(selected_k)}\n")
    return {
        "selected_k": int(selected_k),
        "selection_json": str(json_path),
        "selection_tsv": str(rows_path),
    }

