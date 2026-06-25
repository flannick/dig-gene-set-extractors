from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
import shutil
import subprocess
import sys
from typing import TextIO

from geneset_extractors.preprocessing.rnaseq.cnmf_prepare import (
    _detect_default_threshold,
    _determine_bucket_columns,
    _filter_subset_matrix_for_cnmf,
    _load_meta_index,
    _open_text,
    _parse_csv_list,
    _reservoir_add,
    _resolve_cell_id_index,
    _resolve_gene_id_index,
    _resolve_matrix_orientation,
    _resolve_split_by_cell_type,
    _safe_component,
    _subset_seed,
    _write_meta_subset,
    _write_tmp_subset_from_gene_by_cell,
)
from geneset_extractors.workflows.gtex_runtime_common import write_workflow_provenance_graph


@dataclass
class LigerSubsetPlan:
    subset_id: str
    subset_label: str
    subset_dir: Path
    input_mode: str
    selected_cells: list[str]
    n_cells_before: int
    n_cells_after_downsample: int
    n_cells_after_filter: int = 0
    n_genes_before: int = 0
    n_genes_after_filter: int = 0
    counts_path: Path | None = None
    meta_path: Path | None = None
    run_liger_script_path: Path | None = None
    run_geneset_extractors_from_liger_script_path: Path | None = None
    liger_output_dir: Path | None = None


def _resolve_input_mode(args) -> tuple[str, str]:
    provided: list[tuple[str, str]] = []
    if getattr(args, "matrix_tsv", None):
        provided.append(("matrix_tsv", str(args.matrix_tsv)))
    if getattr(args, "h5ad", None):
        provided.append(("h5ad", str(args.h5ad)))
    if getattr(args, "seurat_rds", None):
        provided.append(("seurat_rds", str(args.seurat_rds)))
    if getattr(args, "mtx_dir", None):
        provided.append(("mtx_dir", str(args.mtx_dir)))
    if len(provided) != 1:
        raise ValueError("Provide exactly one input source: --matrix_tsv, --h5ad, --seurat_rds, or --mtx_dir")
    return provided[0]


def _r_script_path() -> Path:
    return Path(__file__).with_name("liger_inmf.R")


def _write_liger_script(*, subset_dir: Path, input_mode: str, input_path: str, meta_path: str, subset_label: str, args) -> Path:
    script_path = subset_dir / "run_liger.sh"
    out_dir = subset_dir / "liger_out"
    runtime_graph_path = subset_dir / "liger_run.provenance_graph.json"
    prepare_graph_path = subset_dir.parent.parent / "prepare_summary.provenance_graph.json"
    fixed_k = "" if getattr(args, "liger_fixed_k", None) in {None, ""} else str(args.liger_fixed_k)
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        f'R_SCRIPT="{_r_script_path()}"',
        f'INPUT_MODE="{input_mode}"',
        f'INPUT_PATH="{input_path}"',
        f'OUTPUT_DIR="{out_dir}"',
        f'DATASET_COLUMN="{getattr(args, "dataset_column", "") or ""}"',
        f'CELL_TYPE_COLUMN="{getattr(args, "cell_type_column", "") or ""}"',
        f'MAX_CELLS_TOTAL="{int(args.max_cells_total)}"',
        f'MIN_CELLS_PER_CELL_TYPE="{int(args.min_cells_per_cell_type)}"',
        f'SEED="{int(args.seed)}"',
        f'TOP_N_GENES="{int(args.liger_top_n_genes)}"',
        f'META_PATH="{meta_path}"',
        f'CELL_TYPE_LABEL="{subset_label if input_mode == "matrix_tsv" else ""}"',
        f'K_GRID="{str(args.liger_k_grid)}"',
        f'N_REPS="{int(args.liger_n_reps)}"',
        f'FIXED_K="{fixed_k}"',
        f'MIN_CELLS_PER_DATASET="{int(args.liger_min_cells_per_dataset)}"',
        f'MIN_FEATURES="{int(args.liger_min_features)}"',
        f'MIN_UMI="{float(args.liger_min_umi)}"',
        f'MAX_MITO="{float(args.liger_max_mito)}"',
        f'RUNTIME_GRAPH_PATH="{runtime_graph_path}"',
        f'PREPARE_GRAPH_PATH="{prepare_graph_path}"',
        "",
        'Rscript "$R_SCRIPT" "$INPUT_MODE" "$INPUT_PATH" "$OUTPUT_DIR" "$DATASET_COLUMN" "$CELL_TYPE_COLUMN" "$MAX_CELLS_TOTAL" "$MIN_CELLS_PER_CELL_TYPE" "$SEED" "$TOP_N_GENES" "$META_PATH" "$CELL_TYPE_LABEL" "$K_GRID" "$N_REPS" "$FIXED_K" "$MIN_CELLS_PER_DATASET" "$MIN_FEATURES" "$MIN_UMI" "$MAX_MITO"',
        "geneset-extractors workflows scrna_liger_runtime_provenance \\",
        '  --subset_dir "$(pwd)" \\',
        '  --input_mode "$INPUT_MODE" \\',
        '  --input_path "$INPUT_PATH" \\',
        '  --liger_output_dir "$OUTPUT_DIR" \\',
        '  --run_liger_script "$(pwd)/run_liger.sh" \\',
        '  --r_script "$R_SCRIPT" \\',
        '  --runtime_graph_out "$RUNTIME_GRAPH_PATH" \\',
        '  --prepare_provenance_graph_json "$PREPARE_GRAPH_PATH" \\',
        '  --dataset_column "$DATASET_COLUMN" \\',
        '  --cell_type_column "$CELL_TYPE_COLUMN" \\',
        '  --max_cells_total "$MAX_CELLS_TOTAL" \\',
        '  --min_cells_per_cell_type "$MIN_CELLS_PER_CELL_TYPE" \\',
        '  --seed "$SEED" \\',
        '  --liger_top_n_genes "$TOP_N_GENES" \\',
        '  --meta_path "$META_PATH" \\',
        '  --cell_type_label "$CELL_TYPE_LABEL" \\',
        '  --liger_k_grid "$K_GRID" \\',
        '  --liger_n_reps "$N_REPS" \\',
        '  --liger_fixed_k "$FIXED_K" \\',
        '  --liger_min_cells_per_dataset "$MIN_CELLS_PER_DATASET" \\',
        '  --liger_min_features "$MIN_FEATURES" \\',
        '  --liger_min_umi "$MIN_UMI" \\',
        '  --liger_max_mito "$MAX_MITO"',
    ]
    script_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    script_path.chmod(0o755)
    return script_path


def _write_geneset_extractors_from_liger_script(*, subset_dir: Path, args) -> Path:
    script_path = subset_dir / "run_geneset_extractors_from_liger.sh"
    top_k = int(getattr(args, "extractor_top_k", 250))
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        'OUTDIR="$(pwd)/liger_out"',
        'WORKFLOW_GRAPH="$(pwd)/liger_run.provenance_graph.json"',
        "shopt -s nullglob",
        'MATCHES=( "$OUTDIR"/*/gene_loadings.tsv )',
        "shopt -u nullglob",
        'if [[ ${#MATCHES[@]} -eq 0 ]]; then',
        '  echo "error: no LIGER gene_loadings.tsv files found. Run run_liger.sh first." >&2',
        "  exit 2",
        "fi",
        'if [[ ! -f "$WORKFLOW_GRAPH" ]]; then',
        '  echo "error: missing LIGER runtime provenance graph $WORKFLOW_GRAPH. Run run_liger.sh first." >&2',
        "  exit 2",
        "fi",
        'for LOADINGS in "${MATCHES[@]}"; do',
        '  PROGRAM_DIR="$(dirname "$LOADINGS")"',
        '  PROGRAM_LABEL="$(basename "$PROGRAM_DIR")"',
        '  OUT_GENESETS="$PROGRAM_DIR/geneset_extractors_programs"',
        "  geneset-extractors convert rna_sc_programs \\",
        '    --liger_gene_loadings_tsv "$LOADINGS" \\',
        '    --upstream_provenance_graph_json "$WORKFLOW_GRAPH" \\',
        '    --out_dir "$OUT_GENESETS" \\',
        f'    --organism {args.organism} \\',
        f'    --genome_build {args.genome_build} \\',
        "    --score_transform positive \\",
        "    --select top_k \\",
        f"    --top_k {top_k} \\",
        '    --dataset_label "$PROGRAM_LABEL"',
        "done",
    ]
    script_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    script_path.chmod(0o755)
    return script_path


def _execute_script_if_requested(script_path: Path, execute: bool) -> None:
    if not execute:
        return
    if shutil.which("Rscript") is None:
        raise ValueError(
            "Requested --execute true but 'Rscript' was not found on PATH. "
            "Install the required R environment or rerun with --execute false."
        )
    subprocess.run(["bash", str(script_path)], check=True)


def _write_liger_prepare_provenance_graph(
    *,
    args,
    out_dir: Path,
    focus_output_path: Path,
    output_paths: list[tuple[Path, str]],
    input_paths: list[tuple[Path, str]],
) -> Path:
    return write_workflow_provenance_graph(
        workflow_name="scrna_liger_prepare",
        module_name="geneset_extractors.preprocessing.rnaseq.liger_prepare",
        output_dir=out_dir,
        focus_output_path=focus_output_path,
        output_paths=output_paths,
        input_paths=input_paths,
        parameters={
            "input_mode": (
                "matrix_tsv"
                if getattr(args, "matrix_tsv", None)
                else ("h5ad" if getattr(args, "h5ad", None) else ("seurat_rds" if getattr(args, "seurat_rds", None) else "mtx_dir"))
            ),
            "dataset_column": getattr(args, "dataset_column", None),
            "cell_type_column": getattr(args, "cell_type_column", None),
            "split_by_cell_type": getattr(args, "split_by_cell_type", None),
            "max_cells_per_bucket": int(getattr(args, "max_cells_per_bucket", 0) or 0),
            "max_cells_total": int(getattr(args, "max_cells_total", 0) or 0),
            "seed": int(getattr(args, "seed", 0) or 0),
            "liger_k_grid": str(getattr(args, "liger_k_grid", "")),
            "liger_n_reps": int(getattr(args, "liger_n_reps", 0) or 0),
            "liger_fixed_k": getattr(args, "liger_fixed_k", None),
            "liger_top_n_genes": int(getattr(args, "liger_top_n_genes", 0) or 0),
            "liger_min_cells_per_dataset": int(getattr(args, "liger_min_cells_per_dataset", 0) or 0),
            "liger_min_features": int(getattr(args, "liger_min_features", 0) or 0),
            "liger_min_umi": float(getattr(args, "liger_min_umi", 0.0) or 0.0),
            "liger_max_mito": float(getattr(args, "liger_max_mito", 0.0) or 0.0),
            "extractor_top_k": int(getattr(args, "extractor_top_k", 0) or 0),
            "organism": str(getattr(args, "organism", "")),
            "genome_build": str(getattr(args, "genome_build", "")),
        },
        description=(
            "Analysis step that prepares single-cell RNA-seq inputs for LIGER/iNMF, "
            "emits subset manifests and execution scripts, and records the workflow context "
            "used to generate downstream LIGER gene-program loadings."
        ),
    )


def _write_liger_runtime_provenance_graph(
    *,
    subset_dir: Path,
    runtime_graph_out: Path,
    input_mode: str,
    input_path: Path,
    liger_output_dir: Path,
    run_liger_script: Path,
    r_script: Path,
    prepare_provenance_graph_path: Path | None,
    meta_path: Path | None,
    parameters: dict[str, object],
) -> Path:
    output_paths: list[tuple[Path, str]] = []
    for program_dir in sorted(p for p in liger_output_dir.iterdir() if p.is_dir()):
        for filename, role in [
            ("gene_loadings.tsv", "liger_gene_loadings_tsv"),
            ("gene_programs.txt", "liger_gene_programs_txt"),
            ("cell_scores.tsv", "liger_cell_scores_tsv"),
            ("metadata.txt", "liger_metadata_txt"),
            ("k_stability.tsv", "liger_k_stability_tsv"),
            ("factor_importance.txt", "liger_factor_importance_txt"),
        ]:
            candidate = program_dir / filename
            if candidate.exists():
                output_paths.append((candidate, role))
    if not output_paths:
        raise ValueError(f"No LIGER outputs found under {liger_output_dir}")
    focus_output_path = next(
        (path for path, role in output_paths if role == "liger_gene_loadings_tsv"),
        output_paths[0][0],
    )
    input_paths: list[tuple[Path, str]] = [
        (input_path, input_mode),
        (run_liger_script, "run_liger_script"),
        (r_script, "workflow_r_script"),
    ]
    if meta_path is not None and meta_path.exists():
        input_paths.append((meta_path, "meta_tsv"))
    return write_workflow_provenance_graph(
        workflow_name="scrna_liger_runtime_provenance",
        module_name="geneset_extractors.preprocessing.rnaseq.liger_prepare",
        output_dir=subset_dir,
        focus_output_path=focus_output_path,
        output_paths=output_paths,
        input_paths=input_paths,
        parameters=parameters,
        description=(
            "Analysis step that runs LIGER/iNMF on prepared single-cell RNA-seq inputs "
            "and emits gene-program loadings and related latent-factor outputs."
        ),
        upstream_graph_path=prepare_provenance_graph_path,
        graph_path=runtime_graph_out,
    )


def _run_matrix_mode(args, out_dir: Path) -> dict[str, object]:
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

    min_total_per_cell = _detect_default_threshold(getattr(args, "matrix_value_type", "logcounts"), args.min_total_per_cell)
    min_total_per_gene = _detect_default_threshold(getattr(args, "matrix_value_type", "logcounts"), args.min_total_per_gene)

    requested_bucket_columns = _parse_csv_list(args.bucket_columns)
    allow_cell_types = set(_parse_csv_list(args.cell_type_allowlist))

    required_meta_cols = []
    if args.cell_type_column:
        required_meta_cols.append(args.cell_type_column)
    if getattr(args, "dataset_column", None):
        required_meta_cols.append(args.dataset_column)
    required_meta_cols.extend(requested_bucket_columns)
    required_meta_cols = sorted({c for c in required_meta_cols if c})

    meta_fieldnames, meta_by_cell, _n_meta_duplicates = _load_meta_index(
        meta_path=meta_path,
        cell_id_column=args.meta_cell_id_column,
        required_columns=required_meta_cols,
    )

    dataset_active = args.dataset_column if args.dataset_column in meta_fieldnames else None
    cell_type_active = args.cell_type_column if args.cell_type_column in meta_fieldnames else None
    if split_by_cell_type and cell_type_active is None:
        raise ValueError(
            f"split_by_cell_type is enabled but cell_type_column '{args.cell_type_column}' is not present in metadata."
        )

    bucket_columns = _determine_bucket_columns(
        requested_bucket_columns=requested_bucket_columns,
        split_by_cell_type=split_by_cell_type,
        cell_type_column=cell_type_active,
        donor_column=dataset_active,
    )

    import random

    rng = random.Random(int(args.seed))
    bucket_seen: dict[tuple[str, str], int] = {}
    bucket_reservoir: dict[tuple[str, str], list[str]] = {}
    subset_counts_before: dict[str, int] = {}
    n_matrix_rows = 0
    n_matrix_cells_missing_meta = 0
    n_cells_empty_id = 0
    gene_names: list[str] = []
    matrix_cell_id_col_name = ""
    matrix_orientation = "cell_by_gene"
    gene_id_col_idx = 0
    n_genes_input_summary = 0

    with _open_text(matrix_path) as fh:
        reader = csv.reader(fh, delimiter=matrix_delim)
        header = next(reader)
        matrix_orientation = _resolve_matrix_orientation(
            requested=str(args.matrix_orientation),
            header=[str(x) for x in header],
            matrix_cell_id_column=args.matrix_cell_id_column,
            matrix_gene_id_column=args.matrix_gene_id_column,
            meta_cell_ids=set(meta_by_cell.keys()),
        )
        if matrix_orientation == "cell_by_gene":
            cell_idx = _resolve_cell_id_index(header, args.matrix_cell_id_column)
            matrix_cell_id_col_name = header[cell_idx]
            gene_names = [name for i, name in enumerate(header) if i != cell_idx]
            n_genes_input_summary = len(gene_names)
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
        else:
            gene_id_col_idx = _resolve_gene_id_index(header, args.matrix_gene_id_column)
            matrix_cell_id_col_name = args.meta_cell_id_column
            candidate_cells = [str(x).strip() for i, x in enumerate(header) if i != gene_id_col_idx]
            n_matrix_rows = len(candidate_cells)
            for cell_id in candidate_cells:
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
            with _open_text(matrix_path) as fh_count:
                reader_count = csv.reader(fh_count, delimiter=matrix_delim)
                _ = next(reader_count)
                n_genes_input_summary = sum(1 for _ in reader_count)

    retained_subset_ids: list[str] = []
    subset_selected: dict[str, list[str]] = {}
    subset_labels: dict[str, str] = {}
    for subset_id in sorted(subset_counts_before):
        n_before = int(subset_counts_before.get(subset_id, 0))
        if split_by_cell_type and n_before < int(args.min_cells_per_cell_type):
            continue
        cell_ids: list[str] = []
        for (s_id, _bucket_key), reservoir in bucket_reservoir.items():
            if s_id == subset_id:
                cell_ids.extend(reservoir)
        unique_ids = sorted(set(cell_ids))
        if len(unique_ids) > int(args.max_cells_total):
            import random

            subset_rng = random.Random(_subset_seed(int(args.seed), subset_id))
            unique_ids = sorted(subset_rng.sample(unique_ids, int(args.max_cells_total)))
        if unique_ids:
            subset_selected[subset_id] = unique_ids
            retained_subset_ids.append(subset_id)
            subset_labels[subset_id] = subset_id.split("=", 1)[1] if subset_id.startswith("cell_type=") else subset_id

    if not retained_subset_ids:
        raise ValueError("No subsets retained after sampling/filter gates.")

    plans: dict[str, LigerSubsetPlan] = {}
    selected_to_subset: dict[str, str] = {}
    for subset_id in retained_subset_ids:
        subset_dir = subsets_root / subset_id
        subset_dir.mkdir(parents=True, exist_ok=True)
        selected_cells = subset_selected[subset_id]
        for cid in selected_cells:
            selected_to_subset[cid] = subset_id
        plans[subset_id] = LigerSubsetPlan(
            subset_id=subset_id,
            subset_label=subset_labels[subset_id],
            subset_dir=subset_dir,
            input_mode="matrix_tsv",
            selected_cells=selected_cells,
            n_cells_before=int(subset_counts_before.get(subset_id, 0)),
            n_cells_after_downsample=len(selected_cells),
            n_genes_before=int(n_genes_input_summary),
        )

    tmp_paths: dict[str, Path] = {}
    for subset_id, plan in plans.items():
        tmp_path = plan.subset_dir / ".counts_selected.tmp.tsv"
        tmp_paths[subset_id] = tmp_path
        if matrix_orientation == "cell_by_gene":
            with tmp_path.open("w", encoding="utf-8", newline="") as fh:
                writer = csv.writer(fh, delimiter=matrix_delim)
                writer.writerow([matrix_cell_id_col_name] + gene_names)

    if matrix_orientation == "cell_by_gene":
        with _open_text(matrix_path) as fh:
            reader = csv.reader(fh, delimiter=matrix_delim)
            header = next(reader)
            cell_idx = _resolve_cell_id_index(header, args.matrix_cell_id_column)
            gene_indices = [i for i in range(len(header)) if i != cell_idx]
            handle_by_subset: dict[str, TextIO] = {}
            writer_by_subset: dict[str, csv.writer] = {}
            try:
                for subset_id, tmp_path in tmp_paths.items():
                    h = tmp_path.open("a", encoding="utf-8", newline="")
                    handle_by_subset[subset_id] = h
                    writer_by_subset[subset_id] = csv.writer(h, delimiter=matrix_delim)
                for row in reader:
                    if cell_idx >= len(row):
                        continue
                    cell_id = str(row[cell_idx]).strip()
                    subset_id = selected_to_subset.get(cell_id)
                    if not subset_id:
                        continue
                    writer_by_subset[subset_id].writerow([cell_id] + [row[i] if i < len(row) else "" for i in gene_indices])
            finally:
                for h in handle_by_subset.values():
                    h.close()
    else:
        for subset_id in retained_subset_ids:
            plan = plans[subset_id]
            _write_tmp_subset_from_gene_by_cell(
                matrix_path=matrix_path,
                matrix_delim=matrix_delim,
                gene_id_col_idx=gene_id_col_idx,
                selected_cell_ids=plan.selected_cells,
                out_tmp_path=tmp_paths[subset_id],
                out_cell_id_col_name=matrix_cell_id_col_name,
            )

    for subset_id in retained_subset_ids:
        plan = plans[subset_id]
        counts_path = plan.subset_dir / "counts_prefiltered.tsv"
        filter_summary = _filter_subset_matrix_for_cnmf(
            tmp_selected_path=tmp_paths[subset_id],
            out_counts_path=counts_path,
            min_total_per_cell=min_total_per_cell,
            min_total_per_gene=min_total_per_gene,
            matrix_delim=matrix_delim,
            keep_tmp=bool(args.keep_tmp),
        )
        plan.counts_path = counts_path
        plan.n_cells_after_filter = int(filter_summary["n_cells_after_filter"])
        plan.n_genes_after_filter = int(filter_summary["n_genes_after_filter"])
        keep_cells = set(str(x) for x in filter_summary["final_cell_ids"])
        meta_out = plan.subset_dir / "meta.tsv"
        _write_meta_subset(
            meta_path=meta_path,
            out_path=meta_out,
            cell_id_column=args.meta_cell_id_column,
            keep_cells=keep_cells,
        )
        plan.meta_path = meta_out
        plan.liger_output_dir = plan.subset_dir / "liger_out"
        plan.run_liger_script_path = _write_liger_script(
            subset_dir=plan.subset_dir,
            input_mode="matrix_tsv",
            input_path=counts_path.name,
            meta_path=meta_out.name,
            subset_label=plan.subset_label,
            args=args,
        )
        plan.run_geneset_extractors_from_liger_script_path = _write_geneset_extractors_from_liger_script(
            subset_dir=plan.subset_dir,
            args=args,
        )
        _execute_script_if_requested(plan.run_liger_script_path, bool(args.execute))

    manifest_path = out_dir / "subsets_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "subset_id",
                "input_mode",
                "n_cells_before",
                "n_cells_after_downsample",
                "n_cells_after_filter",
                "n_genes_before",
                "n_genes_after_filter",
                "counts_path",
                "meta_path",
                "run_liger_script",
                "run_geneset_extractors_from_liger_script",
            ]
        )
        for subset_id in retained_subset_ids:
            p = plans[subset_id]
            writer.writerow(
                [
                    p.subset_id,
                    p.input_mode,
                    p.n_cells_before,
                    p.n_cells_after_downsample,
                    p.n_cells_after_filter,
                    p.n_genes_before,
                    p.n_genes_after_filter,
                    str((p.counts_path or Path("")).relative_to(out_dir)),
                    str((p.meta_path or Path("")).relative_to(out_dir)),
                    str((p.run_liger_script_path or Path("")).relative_to(out_dir)),
                    str((p.run_geneset_extractors_from_liger_script_path or Path("")).relative_to(out_dir)),
                ]
            )

    summary = {
        "workflow": "scrna_liger_prepare",
        "input_mode": "matrix_tsv",
        "matrix_tsv": str(matrix_path),
        "meta_tsv": str(meta_path),
        "split_by_cell_type": bool(split_by_cell_type),
        "dataset_column": dataset_active,
        "cell_type_column": cell_type_active,
        "bucket_columns": bucket_columns,
        "n_subsets": len(retained_subset_ids),
        "matrix_summary": {
            "n_rows": int(n_matrix_rows),
            "n_genes_input": int(n_genes_input_summary),
            "matrix_orientation": str(matrix_orientation),
            "n_cells_missing_meta": int(n_matrix_cells_missing_meta),
            "n_cells_empty_id": int(n_cells_empty_id),
        },
        "liger": {
            "k_grid": str(args.liger_k_grid),
            "n_reps": int(args.liger_n_reps),
            "fixed_k": getattr(args, "liger_fixed_k", None),
            "top_n_genes": int(args.liger_top_n_genes),
            "execute": bool(args.execute),
        },
        "subsets": [
            {
                "subset_id": plans[sid].subset_id,
                "subset_label": plans[sid].subset_label,
                "counts_path": str((plans[sid].counts_path or Path("")).relative_to(out_dir)),
                "meta_path": str((plans[sid].meta_path or Path("")).relative_to(out_dir)),
                "run_liger_script": str((plans[sid].run_liger_script_path or Path("")).relative_to(out_dir)),
                "run_geneset_extractors_from_liger_script": str(
                    (plans[sid].run_geneset_extractors_from_liger_script_path or Path("")).relative_to(out_dir)
                ),
            }
            for sid in retained_subset_ids
        ],
        "subsets_manifest": str(manifest_path.relative_to(out_dir)),
    }
    summary_path = out_dir / "prepare_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    output_paths: list[tuple[Path, str]] = [
        (summary_path, "workflow_summary"),
        (manifest_path, "subsets_manifest"),
    ]
    for sid in retained_subset_ids:
        plan = plans[sid]
        if plan.counts_path is not None:
            output_paths.append((plan.counts_path, "counts_tsv"))
        if plan.meta_path is not None:
            output_paths.append((plan.meta_path, "meta_tsv"))
        if plan.run_liger_script_path is not None:
            output_paths.append((plan.run_liger_script_path, "run_liger_script"))
        if plan.run_geneset_extractors_from_liger_script_path is not None:
            output_paths.append((plan.run_geneset_extractors_from_liger_script_path, "run_geneset_extractors_script"))
    input_paths: list[tuple[Path, str]] = [
        (matrix_path, "matrix_tsv"),
        (meta_path, "meta_tsv"),
        (_r_script_path(), "workflow_r_script"),
    ]
    graph_path = _write_liger_prepare_provenance_graph(
        args=args,
        out_dir=out_dir,
        focus_output_path=summary_path,
        output_paths=output_paths,
        input_paths=input_paths,
    )
    summary["prepare_provenance_graph_path"] = str(graph_path.relative_to(out_dir))
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(f"prepared scrna_liger subsets={len(retained_subset_ids)} out={out_dir} manifest={manifest_path}", file=sys.stderr)
    return {"out_dir": str(out_dir), "n_subsets": len(retained_subset_ids), "subsets_manifest": str(manifest_path)}


def _run_direct_mode(args, out_dir: Path, input_mode: str, input_path: str) -> dict[str, object]:
    subsets_root = out_dir / "subsets"
    subset_dir = subsets_root / "all"
    subset_dir.mkdir(parents=True, exist_ok=True)
    meta_path = str(getattr(args, "meta_tsv", "") or "")
    run_liger = _write_liger_script(
        subset_dir=subset_dir,
        input_mode=input_mode,
        input_path=input_path,
        meta_path=meta_path,
        subset_label="all",
        args=args,
    )
    run_convert = _write_geneset_extractors_from_liger_script(subset_dir=subset_dir, args=args)
    _execute_script_if_requested(run_liger, bool(args.execute))

    manifest_path = out_dir / "subsets_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "subset_id",
                "input_mode",
                "input_path",
                "meta_path",
                "run_liger_script",
                "run_geneset_extractors_from_liger_script",
            ]
        )
        writer.writerow(
            [
                "all",
                input_mode,
                input_path,
                meta_path,
                str(run_liger.relative_to(out_dir)),
                str(run_convert.relative_to(out_dir)),
            ]
        )
    summary = {
        "workflow": "scrna_liger_prepare",
        "input_mode": input_mode,
        "input_path": input_path,
        "meta_tsv": meta_path or None,
        "dataset_column": getattr(args, "dataset_column", None),
        "cell_type_column": getattr(args, "cell_type_column", None),
        "n_subsets": 1,
        "liger": {
            "k_grid": str(args.liger_k_grid),
            "n_reps": int(args.liger_n_reps),
            "fixed_k": getattr(args, "liger_fixed_k", None),
            "top_n_genes": int(args.liger_top_n_genes),
            "execute": bool(args.execute),
        },
        "subsets": [
            {
                "subset_id": "all",
                "run_liger_script": str(run_liger.relative_to(out_dir)),
                "run_geneset_extractors_from_liger_script": str(run_convert.relative_to(out_dir)),
            }
        ],
        "subsets_manifest": str(manifest_path.relative_to(out_dir)),
    }
    summary_path = out_dir / "prepare_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    output_paths: list[tuple[Path, str]] = [
        (summary_path, "workflow_summary"),
        (manifest_path, "subsets_manifest"),
        (run_liger, "run_liger_script"),
        (run_convert, "run_geneset_extractors_script"),
    ]
    input_paths: list[tuple[Path, str]] = [
        (Path(input_path), input_mode),
        (_r_script_path(), "workflow_r_script"),
    ]
    if meta_path:
        input_paths.append((Path(meta_path), "meta_tsv"))
    graph_path = _write_liger_prepare_provenance_graph(
        args=args,
        out_dir=out_dir,
        focus_output_path=summary_path,
        output_paths=output_paths,
        input_paths=input_paths,
    )
    summary["prepare_provenance_graph_path"] = str(graph_path.relative_to(out_dir))
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(f"prepared scrna_liger subsets=1 out={out_dir} manifest={manifest_path}", file=sys.stderr)
    return {"out_dir": str(out_dir), "n_subsets": 1, "subsets_manifest": str(manifest_path)}


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    input_mode, input_path = _resolve_input_mode(args)
    if input_mode == "matrix_tsv":
        if not getattr(args, "meta_tsv", None):
            raise ValueError("--meta_tsv is required when using --matrix_tsv")
        return _run_matrix_mode(args, out_dir)
    return _run_direct_mode(args, out_dir, input_mode, input_path)


def write_runtime_provenance(args) -> dict[str, object]:
    subset_dir = Path(args.subset_dir).resolve()
    liger_output_dir = Path(args.liger_output_dir)
    if not liger_output_dir.is_absolute():
        liger_output_dir = (subset_dir / liger_output_dir).resolve()
    runtime_graph_out = Path(args.runtime_graph_out)
    if not runtime_graph_out.is_absolute():
        runtime_graph_out = (subset_dir / runtime_graph_out).resolve()
    input_path = Path(args.input_path)
    if not input_path.is_absolute():
        input_path = (subset_dir / input_path).resolve()
    run_liger_script = Path(args.run_liger_script)
    if not run_liger_script.is_absolute():
        run_liger_script = (subset_dir / run_liger_script).resolve()
    r_script = Path(args.r_script)
    if not r_script.is_absolute():
        r_script = (subset_dir / r_script).resolve()
    meta_path: Path | None = None
    raw_meta = str(getattr(args, "meta_path", "") or "").strip()
    if raw_meta:
        candidate = Path(raw_meta)
        meta_path = candidate.resolve() if candidate.is_absolute() else (subset_dir / candidate).resolve()
    prepare_graph: Path | None = None
    raw_prepare_graph = str(getattr(args, "prepare_provenance_graph_json", "") or "").strip()
    if raw_prepare_graph:
        candidate = Path(raw_prepare_graph)
        prepare_graph = candidate.resolve() if candidate.is_absolute() else (subset_dir / candidate).resolve()
    parameters = {
        "input_mode": str(args.input_mode),
        "dataset_column": getattr(args, "dataset_column", None) or None,
        "cell_type_column": getattr(args, "cell_type_column", None) or None,
        "cell_type_label": getattr(args, "cell_type_label", None) or None,
        "max_cells_total": int(args.max_cells_total),
        "min_cells_per_cell_type": int(args.min_cells_per_cell_type),
        "seed": int(args.seed),
        "liger_top_n_genes": int(args.liger_top_n_genes),
        "liger_k_grid": str(args.liger_k_grid),
        "liger_n_reps": int(args.liger_n_reps),
        "liger_fixed_k": getattr(args, "liger_fixed_k", None) or None,
        "liger_min_cells_per_dataset": int(args.liger_min_cells_per_dataset),
        "liger_min_features": int(args.liger_min_features),
        "liger_min_umi": float(args.liger_min_umi),
        "liger_max_mito": float(args.liger_max_mito),
    }
    graph_path = _write_liger_runtime_provenance_graph(
        subset_dir=subset_dir,
        runtime_graph_out=runtime_graph_out,
        input_mode=str(args.input_mode),
        input_path=input_path,
        liger_output_dir=liger_output_dir,
        run_liger_script=run_liger_script,
        r_script=r_script,
        prepare_provenance_graph_path=prepare_graph,
        meta_path=meta_path,
        parameters=parameters,
    )
    return {"runtime_provenance_graph": str(graph_path)}
