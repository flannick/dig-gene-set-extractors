from __future__ import annotations

import csv
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
from typing import TextIO


_SAFE_COMPONENT_RE = re.compile(r"[^A-Za-z0-9._=-]+")


def _parse_float(raw: str) -> float:
    text = str(raw).strip()
    if not text:
        return 0.0
    try:
        value = float(text)
    except ValueError:
        return 0.0
    if value != value:
        return 0.0
    return float(value)


def _open_text(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _safe_component(value: str, fallback: str) -> str:
    out = _SAFE_COMPONENT_RE.sub("_", str(value)).strip("_")
    return out or fallback


def _parse_csv_list(value: str | None) -> list[str]:
    if value is None:
        return []
    out: list[str] = []
    for tok in str(value).split(","):
        token = tok.strip()
        if token:
            out.append(token)
    return out


def _subset_seed(base_seed: int, subset_id: str) -> int:
    digest = hashlib.sha256(f"{base_seed}|{subset_id}".encode("utf-8")).hexdigest()
    return int(digest[:8], 16)


def _reservoir_add(bucket: list[str], seen_count: int, cell_id: str, cap: int, rng: random.Random) -> None:
    if cap <= 0:
        return
    if len(bucket) < cap:
        bucket.append(cell_id)
        return
    j = rng.randint(0, seen_count - 1)
    if j < cap:
        bucket[j] = cell_id


def _detect_default_threshold(value_type: str, user_value: float | None) -> float:
    if user_value is not None:
        return float(user_value)
    if str(value_type) == "counts":
        return 1.0
    return 1e-8


def _resolve_cell_id_index(header: list[str], requested: str | None) -> int:
    if requested:
        req = str(requested).strip()
        if req in header:
            return header.index(req)
        raise ValueError(
            f"matrix_cell_id_column '{req}' was not found in matrix header. "
            f"Available columns: {', '.join(header[:20])}"
        )
    return 0


def _resolve_split_by_cell_type(requested: bool | None, cell_type_column: str | None) -> bool:
    if requested is None:
        return bool(cell_type_column)
    return bool(requested)


@dataclass
class SubsetPlan:
    subset_id: str
    subset_label: str
    subset_dir: Path
    selected_cells: list[str]
    n_cells_before: int
    n_cells_after_downsample: int
    n_cells_written_tmp: int = 0
    n_cells_after_filter: int = 0
    n_cells_dropped_zero_total: int = 0
    n_genes_before: int = 0
    n_genes_after_filter: int = 0
    n_genes_dropped_zero_total: int = 0
    n_non_numeric_values: int = 0
    cnmf_k_list_resolved: list[int] | None = None
    counts_path: Path | None = None
    meta_path: Path | None = None
    run_cnmf_script_path: Path | None = None
    run_cnmf_consensus_auto_k_script_path: Path | None = None
    run_omics2geneset_from_cnmf_script_path: Path | None = None


def _load_meta_index(
    *,
    meta_path: Path,
    cell_id_column: str,
    required_columns: list[str],
) -> tuple[list[str], dict[str, dict[str, str]], int]:
    by_cell: dict[str, dict[str, str]] = {}
    duplicates = 0
    with _open_text(meta_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Metadata table has no header: {meta_path}")
        fieldnames = [str(x) for x in reader.fieldnames]
        if cell_id_column not in fieldnames:
            raise ValueError(
                f"meta_cell_id_column '{cell_id_column}' not found in metadata header. "
                f"Available columns: {', '.join(fieldnames)}"
            )
        for col in required_columns:
            if col and col not in fieldnames:
                raise ValueError(
                    f"Required metadata column '{col}' not found in metadata header. "
                    f"Available columns: {', '.join(fieldnames)}"
                )
        for row in reader:
            cell_id = str(row.get(cell_id_column, "")).strip()
            if not cell_id:
                continue
            if cell_id in by_cell:
                duplicates += 1
            keep = {cell_id_column: cell_id}
            for col in required_columns:
                if not col:
                    continue
                keep[col] = str(row.get(col, "")).strip()
            by_cell[cell_id] = keep
    return fieldnames, by_cell, duplicates


def _determine_bucket_columns(
    *,
    requested_bucket_columns: list[str],
    split_by_cell_type: bool,
    cell_type_column: str | None,
    donor_column: str | None,
) -> list[str]:
    if requested_bucket_columns:
        return requested_bucket_columns
    cols: list[str] = []
    if split_by_cell_type and cell_type_column:
        cols.append(cell_type_column)
        if donor_column:
            cols.append(donor_column)
        return cols
    if donor_column:
        cols.append(donor_column)
    return cols


def _write_meta_subset(
    *,
    meta_path: Path,
    out_path: Path,
    cell_id_column: str,
    keep_cells: set[str],
) -> None:
    with _open_text(meta_path) as in_fh:
        reader = csv.DictReader(in_fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Metadata table has no header: {meta_path}")
        fieldnames = [str(x) for x in reader.fieldnames]
        with out_path.open("w", encoding="utf-8", newline="") as out_fh:
            writer = csv.DictWriter(out_fh, delimiter="\t", fieldnames=fieldnames)
            writer.writeheader()
            for row in reader:
                cell_id = str(row.get(cell_id_column, "")).strip()
                if cell_id and cell_id in keep_cells:
                    writer.writerow({k: ("" if v is None else str(v)) for k, v in row.items()})


def _filter_subset_matrix_for_cnmf(
    *,
    tmp_selected_path: Path,
    out_counts_path: Path,
    min_total_per_cell: float,
    min_total_per_gene: float,
    matrix_delim: str,
    keep_tmp: bool,
) -> dict[str, object]:
    with tmp_selected_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter=matrix_delim)
        header = next(reader)
    if len(header) < 2:
        raise ValueError(f"Subset temp matrix has no gene columns: {tmp_selected_path}")

    current_input = tmp_selected_path
    current_genes = list(header[1:])
    n_genes_before = len(current_genes)
    total_cells_dropped = 0
    total_genes_dropped = 0
    n_non_numeric = 0
    keep_cell_ids: list[str] = []
    stable_output: Path | None = None

    temp_paths: list[Path] = []
    max_iter = 6
    for iter_idx in range(1, max_iter + 1):
        iter_out = out_counts_path.parent / f".counts_iter_{iter_idx}.tsv"
        temp_paths.append(iter_out)
        gene_sums = [0.0 for _ in current_genes]
        n_cells_in = 0
        n_cells_out = 0
        keep_cell_ids_iter: list[str] = []

        with current_input.open("r", encoding="utf-8", newline="") as in_fh, iter_out.open(
            "w", encoding="utf-8", newline=""
        ) as out_fh:
            reader = csv.reader(in_fh, delimiter=matrix_delim)
            header_in = next(reader)
            cell_col = header_in[0]
            gene_index = {name: i for i, name in enumerate(header_in)}
            target_indices = [gene_index[g] for g in current_genes if g in gene_index]
            writer = csv.writer(out_fh, delimiter=matrix_delim)
            writer.writerow([cell_col] + current_genes)

            for row in reader:
                n_cells_in += 1
                if not row:
                    continue
                cell_id = str(row[0]).strip()
                values: list[float] = []
                cell_total = 0.0
                for gi, col_idx in enumerate(target_indices):
                    raw = row[col_idx] if col_idx < len(row) else ""
                    val = _parse_float(raw)
                    if str(raw).strip() and val == 0.0 and str(raw).strip() not in {"0", "0.0", "0.00"}:
                        n_non_numeric += 1
                    values.append(val)
                    cell_total += val
                if cell_total < float(min_total_per_cell):
                    continue
                n_cells_out += 1
                keep_cell_ids_iter.append(cell_id)
                writer.writerow([cell_id] + values)
                for gi, val in enumerate(values):
                    gene_sums[gi] += val

        total_cells_dropped += max(n_cells_in - n_cells_out, 0)
        new_genes = [
            g
            for g, total in zip(current_genes, gene_sums)
            if float(total) >= float(min_total_per_gene) and float(total) > 0.0
        ]
        total_genes_dropped += max(len(current_genes) - len(new_genes), 0)
        if not new_genes:
            raise ValueError(
                "No genes remained after cNMF prefiltering. "
                "Lower --min_total_per_gene and/or increase downsampled cells."
            )

        if new_genes == current_genes:
            stable_output = iter_out
            keep_cell_ids = keep_cell_ids_iter
            break

        current_input = iter_out
        current_genes = new_genes

    if stable_output is None:
        raise ValueError(
            "cNMF prefiltering did not converge within iteration limit. "
            "Check input matrix and threshold settings."
        )

    if not keep_cell_ids:
        raise ValueError(
            "No cells remained after cNMF prefiltering. "
            "Lower --min_total_per_cell and/or increase downsampled cells."
        )

    stable_output.replace(out_counts_path)

    if not keep_tmp and tmp_selected_path.exists():
        tmp_selected_path.unlink()
    if not keep_tmp:
        for p in temp_paths:
            if p.exists() and p != out_counts_path:
                p.unlink()

    return {
        "n_cells_after_filter": len(keep_cell_ids),
        "n_cells_dropped_zero_total": total_cells_dropped,
        "n_genes_before": n_genes_before,
        "n_genes_after_filter": len(current_genes),
        "n_genes_dropped_zero_total": total_genes_dropped,
        "n_non_numeric_values": n_non_numeric,
        "final_gene_names": current_genes,
        "final_cell_ids": keep_cell_ids,
    }


def _parse_k_list(value: str) -> list[int]:
    tokens = [tok for tok in str(value).replace(",", " ").split() if tok]
    if not tokens:
        raise ValueError("cnmf_k_list must contain at least one K value")
    out: list[int] = []
    for token in tokens:
        try:
            k = int(token)
        except ValueError as exc:
            raise ValueError(f"Invalid K token in cnmf_k_list: {token}") from exc
        if k <= 0:
            raise ValueError(f"K values must be > 0, got: {k}")
        if k not in out:
            out.append(k)
    return out


def _resolve_cnmf_k_list_tokens(cnmf_k_list: str, n_cells_after_downsample: int) -> tuple[list[int], str]:
    raw = str(cnmf_k_list).strip().lower()
    if raw != "auto":
        return _parse_k_list(cnmf_k_list), "explicit"
    n = int(n_cells_after_downsample)
    if n <= 1500:
        return [5, 8, 10, 12, 15], "auto_small"
    if n <= 5000:
        return [10, 15, 20, 25, 30], "auto_medium"
    return [15, 20, 25, 30, 35, 40, 45, 50], "auto_large"


def _parse_fixed_k_or_auto(raw_value: str) -> tuple[str, int | None]:
    text = str(raw_value).strip().lower()
    if text == "auto":
        return "auto", None
    try:
        k = int(text)
    except ValueError as exc:
        raise ValueError(f"--cnmf_k must be 'auto' or an integer, got: {raw_value}") from exc
    if k <= 0:
        raise ValueError(f"--cnmf_k must be > 0 when integer, got: {k}")
    return "fixed", k


def _selector_call_line(args, quote_fixed_k: str | None = None) -> str:
    cmd = (
        'omics2geneset workflows cnmf_select_k --cnmf_output_dir "$OUTDIR" --name "$NAME" '
        f'--strategy {args.cnmf_select_strategy} '
        f'--stability_frac_of_max {float(args.cnmf_select_stability_frac_of_max)} '
        f'--min_stability_abs {float(args.cnmf_select_min_stability_abs)} '
        f'--require_local_max {"true" if bool(args.cnmf_select_require_local_max) else "false"}'
    )
    if quote_fixed_k is not None:
        cmd += f" --fixed_k {quote_fixed_k}"
    return cmd


def _selector_manual_call_line(quote_fixed_k: str) -> str:
    return (
        'omics2geneset workflows cnmf_select_k --cnmf_output_dir "$OUTDIR" --name "$NAME" '
        f"--strategy manual --fixed_k {quote_fixed_k}"
    )


def _write_cnmf_prepare_script(
    *,
    subset_dir: Path,
    subset_id: str,
    args,
    counts_filename: str,
    k_tokens: list[int],
) -> Path:
    script_path = subset_dir / "run_cnmf.sh"
    safe_subset = _safe_component(subset_id, "subset")
    cnmf_name = str(args.cnmf_name).strip() if args.cnmf_name else safe_subset
    densify = bool(args.cnmf_densify)

    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        'OUTDIR="$(pwd)/cnmf_out"',
        f'NAME="{cnmf_name}"',
        f'COUNTS="{counts_filename}"',
        f'K_LIST="{" ".join(str(k) for k in k_tokens)}"',
        f'N_ITER="{int(args.cnmf_n_iter)}"',
        f'NUMGENES="{int(args.cnmf_numgenes)}"',
        f'SEED="{int(args.seed)}"',
        f'TOTAL_WORKERS="{int(args.cnmf_total_workers)}"',
        "",
        'mkdir -p "$OUTDIR"',
        "",
        "cnmf prepare \\",
        '  --output-dir "$OUTDIR" \\',
        '  --name "$NAME" \\',
        '  --counts "$COUNTS" \\',
        '  --components $K_LIST \\',
        '  --n-iter "$N_ITER" \\',
        '  --numgenes "$NUMGENES" \\',
        '  --seed "$SEED" \\',
        ("  --densify" if densify else "  # --densify  # enable if your cNMF install expects dense text"),
        "",
        "cnmf factorize \\",
        '  --output-dir "$OUTDIR" \\',
        '  --name "$NAME" \\',
        '  --worker-index 0 \\',
        '  --total-workers "$TOTAL_WORKERS"',
        "",
        'cnmf combine --output-dir "$OUTDIR" --name "$NAME"',
        'cnmf k_selection_plot --output-dir "$OUTDIR" --name "$NAME"',
    ]
    script_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    script_path.chmod(0o755)
    return script_path


def _write_cnmf_consensus_auto_k_script(
    *,
    subset_dir: Path,
    subset_id: str,
    args,
) -> Path:
    script_path = subset_dir / "run_cnmf_consensus_auto_k.sh"
    safe_subset = _safe_component(subset_id, "subset")
    cnmf_name = str(args.cnmf_name).strip() if args.cnmf_name else safe_subset
    _k_mode, fixed_k = _parse_fixed_k_or_auto(args.cnmf_k)
    selector_cmd_auto = _selector_call_line(args)
    selector_cmd_manual = _selector_manual_call_line('"$K"')

    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        'OUTDIR="$(pwd)/cnmf_out"',
        f'NAME="{cnmf_name}"',
        f'K_OVERRIDE="{args.cnmf_k}"',
        "",
        'if [[ "$K_OVERRIDE" == "auto" ]]; then',
        f'  K="$({selector_cmd_auto})"',
        "else",
        '  K="$K_OVERRIDE"',
        (f"  {selector_cmd_manual} > /dev/null" if fixed_k is not None else f"  {selector_cmd_auto} > /dev/null"),
        "fi",
        "",
        'echo "Using K=$K for consensus" >&2',
        "cnmf consensus \\",
        '  --output-dir "$OUTDIR" \\',
        '  --name "$NAME" \\',
        '  --components "$K" \\',
        f'  --local-density-threshold {float(args.cnmf_local_density_threshold)} \\',
        f'  --local-neighborhood-size {float(args.cnmf_local_neighborhood_size)} \\',
        ("  --show-clustering" if bool(args.cnmf_show_clustering) else "  # --show-clustering"),
    ]
    script_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    script_path.chmod(0o755)
    return script_path


def _write_omics2geneset_from_cnmf_script(
    *,
    subset_dir: Path,
    subset_id: str,
    args,
) -> Path:
    script_path = subset_dir / "run_omics2geneset_from_cnmf.sh"
    safe_subset = _safe_component(subset_id, "subset")
    cnmf_name = str(args.cnmf_name).strip() if args.cnmf_name else safe_subset
    export_kind = str(args.cnmf_export_kind).strip()
    if export_kind not in {"tpm", "score"}:
        raise ValueError(f"Unsupported --cnmf_export_kind: {export_kind}")

    extra_flags = ""
    if export_kind == "score":
        extra_flags = "--score_transform signed --gmt_split_signed true"
    selector_cmd_auto = _selector_call_line(args)
    selector_cmd_manual = _selector_manual_call_line('"$K"')
    k_is_auto = str(args.cnmf_k).strip().lower() == "auto"

    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        'OUTDIR="$(pwd)/cnmf_out"',
        f'NAME="{cnmf_name}"',
        f'K_OVERRIDE="{args.cnmf_k}"',
        f'EXPORT_KIND="{export_kind}"',
        "",
        'if [[ "$K_OVERRIDE" == "auto" ]]; then',
        f'  K="$({selector_cmd_auto})"',
        "else",
        '  K="$K_OVERRIDE"',
        (f"  {selector_cmd_manual} > /dev/null" if not k_is_auto else f"  {selector_cmd_auto} > /dev/null"),
        "fi",
        "",
        "shopt -s nullglob",
        'MATCHES=( "$OUTDIR/$NAME.gene_spectra_${EXPORT_KIND}.k_${K}.dt_"*.txt )',
        "shopt -u nullglob",
        'if [[ ${#MATCHES[@]} -eq 0 ]]; then',
        '  echo "error: no gene spectra file matched kind=${EXPORT_KIND}, k=${K}. Run consensus first." >&2',
        "  exit 2",
        "fi",
        'if [[ ${#MATCHES[@]} -gt 1 ]]; then',
        '  echo "error: multiple spectra files matched for kind=${EXPORT_KIND}, k=${K}:" >&2',
        '  printf "  %s\\n" "${MATCHES[@]}" >&2',
        '  echo "Please keep one file or edit this script to pick a specific dt." >&2',
        "  exit 2",
        "fi",
        'SPECTRA="${MATCHES[0]}"',
        'OUT_GENESETS="$OUTDIR/omics2geneset_programs_k_${K}_${EXPORT_KIND}"',
        "",
        "omics2geneset convert rna_sc_programs \\",
        '  --cnmf_gene_spectra_tsv "$SPECTRA" \\',
        '  --out_dir "$OUT_GENESETS" \\',
        f'  --organism {args.organism} --genome_build {args.genome_build}' + (f" {extra_flags}" if extra_flags else ""),
    ]
    script_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    script_path.chmod(0o755)
    return script_path


def _execute_script_if_requested(script_path: Path, execute: bool) -> None:
    if not execute:
        return
    if shutil.which("cnmf") is None:
        raise ValueError(
            "Requested --execute true but 'cnmf' was not found on PATH. "
            "Install cNMF or rerun with --execute false."
        )
    subprocess.run(["bash", str(script_path)], check=True)


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    subsets_root = out_dir / "subsets"
    subsets_root.mkdir(parents=True, exist_ok=True)

    matrix_path = Path(args.matrix_tsv)
    meta_path = Path(args.meta_tsv)
    matrix_delim = str(args.matrix_delim)
    if len(matrix_delim) != 1:
        raise ValueError("--matrix_delim must be a single character delimiter")

    split_by_cell_type = _resolve_split_by_cell_type(args.split_by_cell_type, args.cell_type_column)
    if split_by_cell_type and not args.cell_type_column:
        raise ValueError("split_by_cell_type requires --cell_type_column")

    min_total_per_cell = _detect_default_threshold(args.matrix_value_type, args.min_total_per_cell)
    min_total_per_gene = _detect_default_threshold(args.matrix_value_type, args.min_total_per_gene)
    if str(args.cnmf_k).strip().lower() == "auto" and str(args.cnmf_select_strategy).strip().lower() == "manual":
        raise ValueError(
            "--cnmf_select_strategy manual requires explicit --cnmf_k <INT> for generated scripts "
            "(or use a non-manual strategy with --cnmf_k auto)."
        )

    requested_bucket_columns = _parse_csv_list(args.bucket_columns)
    allow_cell_types = set(_parse_csv_list(args.cell_type_allowlist))

    required_meta_cols = []
    if args.cell_type_column:
        required_meta_cols.append(args.cell_type_column)
    if args.donor_column:
        required_meta_cols.append(args.donor_column)
    if args.phenotype_column:
        required_meta_cols.append(args.phenotype_column)
    required_meta_cols.extend(requested_bucket_columns)
    required_meta_cols = sorted({c for c in required_meta_cols if c})

    meta_fieldnames, meta_by_cell, n_meta_duplicates = _load_meta_index(
        meta_path=meta_path,
        cell_id_column=args.meta_cell_id_column,
        required_columns=required_meta_cols,
    )
    if n_meta_duplicates > 0:
        print(
            f"warning: metadata had duplicate cell_id rows; last occurrence used (n_duplicates={n_meta_duplicates}).",
            file=sys.stderr,
        )

    donor_active = args.donor_column if args.donor_column in meta_fieldnames else None
    if args.donor_column and donor_active is None:
        print(
            f"warning: donor_column '{args.donor_column}' not found in metadata; donor-aware bucketing disabled.",
            file=sys.stderr,
        )
    cell_type_active = args.cell_type_column if args.cell_type_column in meta_fieldnames else None
    if split_by_cell_type and cell_type_active is None:
        raise ValueError(
            f"split_by_cell_type is enabled but cell_type_column '{args.cell_type_column}' is not present in metadata."
        )

    bucket_columns = _determine_bucket_columns(
        requested_bucket_columns=requested_bucket_columns,
        split_by_cell_type=split_by_cell_type,
        cell_type_column=cell_type_active,
        donor_column=donor_active,
    )

    rng = random.Random(int(args.seed))
    bucket_seen: dict[tuple[str, str], int] = {}
    bucket_reservoir: dict[tuple[str, str], list[str]] = {}
    subset_counts_before: dict[str, int] = {}

    n_matrix_rows = 0
    n_matrix_cells_missing_meta = 0
    n_matrix_cells_skipped_allowlist = 0
    n_cells_empty_id = 0
    gene_names: list[str] = []
    matrix_cell_id_col_name = ""

    with _open_text(matrix_path) as fh:
        reader = csv.reader(fh, delimiter=matrix_delim)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"matrix_tsv is empty: {matrix_path}") from exc
        if not header:
            raise ValueError(f"matrix_tsv has empty header: {matrix_path}")

        cell_idx = _resolve_cell_id_index(header, args.matrix_cell_id_column)
        matrix_cell_id_col_name = header[cell_idx]
        gene_names = [name for i, name in enumerate(header) if i != cell_idx]
        if not gene_names:
            raise ValueError("matrix_tsv must include at least one gene column")

        for row in reader:
            n_matrix_rows += 1
            if cell_idx >= len(row):
                continue
            cell_id = str(row[cell_idx]).strip()
            if not cell_id:
                n_cells_empty_id += 1
                continue
            meta = meta_by_cell.get(cell_id)
            if meta is None:
                n_matrix_cells_missing_meta += 1
                continue

            if split_by_cell_type:
                raw_cell_type = str(meta.get(cell_type_active or "", "")).strip() if cell_type_active else ""
                subset_label = raw_cell_type or "unknown"
                if allow_cell_types and subset_label not in allow_cell_types:
                    n_matrix_cells_skipped_allowlist += 1
                    continue
                subset_id = f"cell_type={_safe_component(subset_label, 'unknown')}"
            else:
                subset_label = "all"
                subset_id = "all"

            subset_counts_before[subset_id] = int(subset_counts_before.get(subset_id, 0)) + 1

            if bucket_columns:
                bucket_vals = [str(meta.get(col, "")).strip() or "unknown" for col in bucket_columns]
                bucket_key = "|".join(bucket_vals)
            else:
                bucket_key = "__all__"

            skey = (subset_id, bucket_key)
            seen = int(bucket_seen.get(skey, 0)) + 1
            bucket_seen[skey] = seen
            bucket = bucket_reservoir.setdefault(skey, [])
            _reservoir_add(bucket, seen, cell_id, int(args.max_cells_per_bucket), rng)

    retained_subset_ids: list[str] = []
    subset_selected: dict[str, list[str]] = {}
    subset_labels: dict[str, str] = {}
    for subset_id in sorted(subset_counts_before):
        n_before = int(subset_counts_before.get(subset_id, 0))
        if split_by_cell_type and n_before < int(args.min_cells_per_cell_type):
            print(
                "warning: skipping subset due to min_cells_per_cell_type gate: "
                f"{subset_id} n_cells_before={n_before} < {int(args.min_cells_per_cell_type)}",
                file=sys.stderr,
            )
            continue

        cell_ids: list[str] = []
        for (s_id, _bucket_key), reservoir in bucket_reservoir.items():
            if s_id != subset_id:
                continue
            cell_ids.extend(reservoir)
        unique_ids = sorted(set(cell_ids))
        if not unique_ids:
            print(f"warning: subset {subset_id} has no selected cells after bucket sampling; skipping.", file=sys.stderr)
            continue

        if len(unique_ids) > int(args.max_cells_total):
            subset_rng = random.Random(_subset_seed(int(args.seed), subset_id))
            unique_ids = sorted(subset_rng.sample(unique_ids, int(args.max_cells_total)))

        subset_selected[subset_id] = unique_ids
        retained_subset_ids.append(subset_id)
        subset_labels[subset_id] = subset_id.split("=", 1)[1] if subset_id.startswith("cell_type=") else subset_id

    if not retained_subset_ids:
        raise ValueError(
            "No subsets retained after sampling/filter gates. "
            "Adjust min_cells_per_cell_type, allowlist, or sampling caps."
        )

    plans: dict[str, SubsetPlan] = {}
    selected_to_subset: dict[str, str] = {}
    for subset_id in retained_subset_ids:
        subset_dir = subsets_root / subset_id
        subset_dir.mkdir(parents=True, exist_ok=True)
        selected_cells = subset_selected[subset_id]
        for cid in selected_cells:
            selected_to_subset[cid] = subset_id
        plans[subset_id] = SubsetPlan(
            subset_id=subset_id,
            subset_label=subset_labels[subset_id],
            subset_dir=subset_dir,
            selected_cells=selected_cells,
            n_cells_before=int(subset_counts_before.get(subset_id, 0)),
            n_cells_after_downsample=len(selected_cells),
            n_genes_before=len(gene_names),
        )

    tmp_handles: dict[str, TextIO] = {}
    tmp_writers: dict[str, csv.writer] = {}
    tmp_paths: dict[str, Path] = {}
    for subset_id, plan in plans.items():
        tmp_path = plan.subset_dir / ".counts_selected.tmp.tsv"
        tmp_paths[subset_id] = tmp_path
        fh = tmp_path.open("w", encoding="utf-8", newline="")
        writer = csv.writer(fh, delimiter=matrix_delim)
        writer.writerow([matrix_cell_id_col_name] + gene_names)
        tmp_handles[subset_id] = fh
        tmp_writers[subset_id] = writer

    with _open_text(matrix_path) as fh:
        reader = csv.reader(fh, delimiter=matrix_delim)
        header = next(reader)
        cell_idx = _resolve_cell_id_index(header, args.matrix_cell_id_column)
        gene_indices = [i for i in range(len(header)) if i != cell_idx]
        for row in reader:
            if cell_idx >= len(row):
                continue
            cell_id = str(row[cell_idx]).strip()
            subset_id = selected_to_subset.get(cell_id)
            if not subset_id:
                continue
            out_row = [cell_id] + [row[i] if i < len(row) else "" for i in gene_indices]
            tmp_writers[subset_id].writerow(out_row)
            plans[subset_id].n_cells_written_tmp += 1

    for fh in tmp_handles.values():
        fh.close()

    for subset_id in retained_subset_ids:
        plan = plans[subset_id]
        tmp_path = tmp_paths[subset_id]
        counts_path = plan.subset_dir / "counts_prefiltered.tsv"
        filter_summary = _filter_subset_matrix_for_cnmf(
            tmp_selected_path=tmp_path,
            out_counts_path=counts_path,
            min_total_per_cell=min_total_per_cell,
            min_total_per_gene=min_total_per_gene,
            matrix_delim=matrix_delim,
            keep_tmp=bool(args.keep_tmp),
        )
        plan.counts_path = counts_path
        plan.n_cells_after_filter = int(filter_summary["n_cells_after_filter"])
        plan.n_cells_dropped_zero_total = int(filter_summary["n_cells_dropped_zero_total"])
        plan.n_genes_after_filter = int(filter_summary["n_genes_after_filter"])
        plan.n_genes_dropped_zero_total = int(filter_summary["n_genes_dropped_zero_total"])
        plan.n_non_numeric_values = int(filter_summary["n_non_numeric_values"])

        keep_cells = set(str(x) for x in filter_summary["final_cell_ids"])
        meta_out = plan.subset_dir / "meta.tsv"
        _write_meta_subset(
            meta_path=meta_path,
            out_path=meta_out,
            cell_id_column=args.meta_cell_id_column,
            keep_cells=keep_cells,
        )
        plan.meta_path = meta_out

        if plan.n_cells_dropped_zero_total > 0:
            print(
                f"warning: dropped {plan.n_cells_dropped_zero_total} zero-total cells in subset {subset_id}. "
                "Adjust --min_total_per_cell if this is unexpected.",
                file=sys.stderr,
            )
        if plan.n_genes_dropped_zero_total > 0:
            print(
                f"warning: dropped {plan.n_genes_dropped_zero_total} zero-total/low-total genes in subset {subset_id}. "
                "Adjust --min_total_per_gene if this is unexpected.",
                file=sys.stderr,
            )

        k_tokens, _k_mode = _resolve_cnmf_k_list_tokens(args.cnmf_k_list, plan.n_cells_after_downsample)
        plan.cnmf_k_list_resolved = k_tokens
        plan.run_cnmf_script_path = _write_cnmf_prepare_script(
            subset_dir=plan.subset_dir,
            subset_id=subset_id,
            args=args,
            counts_filename=counts_path.name,
            k_tokens=k_tokens,
        )
        plan.run_cnmf_consensus_auto_k_script_path = _write_cnmf_consensus_auto_k_script(
            subset_dir=plan.subset_dir,
            subset_id=subset_id,
            args=args,
        )
        plan.run_omics2geneset_from_cnmf_script_path = _write_omics2geneset_from_cnmf_script(
            subset_dir=plan.subset_dir,
            subset_id=subset_id,
            args=args,
        )
        _execute_script_if_requested(plan.run_cnmf_script_path, bool(args.execute))

    manifest_path = out_dir / "subsets_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "subset_id",
                "n_cells_before",
                "n_cells_after_downsample",
                "n_cells_after_filter",
                "n_genes_before",
                "n_genes_after_filter",
                "cnmf_k_list_resolved",
                "counts_path",
                "meta_path",
                "run_cnmf_script",
                "run_cnmf_consensus_auto_k_script",
                "run_omics2geneset_from_cnmf_script",
            ]
        )
        for subset_id in retained_subset_ids:
            p = plans[subset_id]
            writer.writerow(
                [
                    p.subset_id,
                    p.n_cells_before,
                    p.n_cells_after_downsample,
                    p.n_cells_after_filter,
                    p.n_genes_before,
                    p.n_genes_after_filter,
                    " ".join(str(x) for x in (p.cnmf_k_list_resolved or [])),
                    str((p.counts_path or Path("")).relative_to(out_dir)),
                    str((p.meta_path or Path("")).relative_to(out_dir)),
                    str((p.run_cnmf_script_path or Path("")).relative_to(out_dir)),
                    str((p.run_cnmf_consensus_auto_k_script_path or Path("")).relative_to(out_dir)),
                    str((p.run_omics2geneset_from_cnmf_script_path or Path("")).relative_to(out_dir)),
                ]
            )

    summary = {
        "workflow": "scrna_cnmf_prepare",
        "matrix_tsv": str(matrix_path),
        "meta_tsv": str(meta_path),
        "matrix_value_type": str(args.matrix_value_type),
        "seed": int(args.seed),
        "split_by_cell_type": bool(split_by_cell_type),
        "cell_type_column": args.cell_type_column,
        "donor_column": donor_active,
        "phenotype_column": args.phenotype_column if args.phenotype_column in meta_fieldnames else None,
        "bucket_columns": bucket_columns,
        "cell_type_allowlist": sorted(allow_cell_types) if allow_cell_types else [],
        "min_cells_per_cell_type": int(args.min_cells_per_cell_type),
        "max_cells_per_bucket": int(args.max_cells_per_bucket),
        "max_cells_total": int(args.max_cells_total),
        "min_total_per_cell": float(min_total_per_cell),
        "min_total_per_gene": float(min_total_per_gene),
        "cnmf": {
            "cnmf_name": args.cnmf_name,
            "cnmf_k_list_input": str(args.cnmf_k_list),
            "cnmf_k_input": str(args.cnmf_k),
            "cnmf_select_strategy": str(args.cnmf_select_strategy),
            "cnmf_select_stability_frac_of_max": float(args.cnmf_select_stability_frac_of_max),
            "cnmf_select_min_stability_abs": float(args.cnmf_select_min_stability_abs),
            "cnmf_select_require_local_max": bool(args.cnmf_select_require_local_max),
            "cnmf_n_iter": int(args.cnmf_n_iter),
            "cnmf_numgenes": int(args.cnmf_numgenes),
            "cnmf_total_workers": int(args.cnmf_total_workers),
            "cnmf_densify": bool(args.cnmf_densify),
            "cnmf_local_density_threshold": float(args.cnmf_local_density_threshold),
            "cnmf_local_neighborhood_size": float(args.cnmf_local_neighborhood_size),
            "cnmf_show_clustering": bool(args.cnmf_show_clustering),
            "cnmf_export_kind": str(args.cnmf_export_kind),
            "execute": bool(args.execute),
        },
        "matrix_summary": {
            "n_rows": int(n_matrix_rows),
            "n_genes_input": int(len(gene_names)),
            "n_cells_missing_meta": int(n_matrix_cells_missing_meta),
            "n_cells_empty_id": int(n_cells_empty_id),
            "n_cells_skipped_allowlist": int(n_matrix_cells_skipped_allowlist),
        },
        "n_subsets": len(retained_subset_ids),
        "subsets": [
            {
                "subset_id": plans[sid].subset_id,
                "subset_label": plans[sid].subset_label,
                "n_cells_before": plans[sid].n_cells_before,
                "n_cells_after_downsample": plans[sid].n_cells_after_downsample,
                "n_cells_after_filter": plans[sid].n_cells_after_filter,
                "n_cells_dropped_zero_total": plans[sid].n_cells_dropped_zero_total,
                "n_genes_before": plans[sid].n_genes_before,
                "n_genes_after_filter": plans[sid].n_genes_after_filter,
                "n_genes_dropped_zero_total": plans[sid].n_genes_dropped_zero_total,
                "n_non_numeric_values": plans[sid].n_non_numeric_values,
                "cnmf_k_list_resolved": plans[sid].cnmf_k_list_resolved or [],
                "counts_path": str((plans[sid].counts_path or Path("")).relative_to(out_dir)),
                "meta_path": str((plans[sid].meta_path or Path("")).relative_to(out_dir)),
                "run_cnmf_script": str((plans[sid].run_cnmf_script_path or Path("")).relative_to(out_dir)),
                "run_cnmf_consensus_auto_k_script": str(
                    (plans[sid].run_cnmf_consensus_auto_k_script_path or Path("")).relative_to(out_dir)
                ),
                "run_omics2geneset_from_cnmf_script": str(
                    (plans[sid].run_omics2geneset_from_cnmf_script_path or Path("")).relative_to(out_dir)
                ),
            }
            for sid in retained_subset_ids
        ],
        "subsets_manifest": str(manifest_path.relative_to(out_dir)),
    }
    (out_dir / "prepare_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    print(
        "prepared scrna_cnmf subsets="
        f"{len(retained_subset_ids)} out={out_dir} manifest={manifest_path}",
        file=sys.stderr,
    )
    return {
        "out_dir": str(out_dir),
        "n_subsets": len(retained_subset_ids),
        "subsets_manifest": str(manifest_path),
    }
