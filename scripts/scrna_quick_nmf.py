#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _load_optional_deps():
    try:
        import numpy as np  # type: ignore
        import pandas as pd  # type: ignore
        from sklearn.decomposition import NMF  # type: ignore
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised by user env, not unit tests
        raise SystemExit(
            "Missing optional dependency for scrna_quick_nmf.py: "
            f"{exc}. Install optional extras with: python -m pip install -e '.[scrna_tools]'"
        )
    return np, pd, NMF


def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Quick helper: downsample a scRNA matrix and fit sklearn NMF, "
            "then write genes-by-program loadings TSV for omics2geneset rna_sc_programs."
        )
    )
    ap.add_argument("--matrix_tsv", required=True, help="Cell x gene matrix TSV (rows=cells, columns=genes)")
    ap.add_argument("--out_tsv", required=True, help="Output genes-by-program loadings TSV")
    ap.add_argument("--metadata_tsv", help="Optional metadata TSV with at least cell_id")
    ap.add_argument("--cell_id_column", default="cell_id", help="Cell identifier column name")
    ap.add_argument("--cell_type_column", default="cell_type", help="Optional cell-type column for bucketed downsampling")
    ap.add_argument("--donor_id_column", default="donor_id", help="Optional donor column for bucketed downsampling")
    ap.add_argument("--max_cells_total", type=int, default=20000, help="Maximum total cells after downsampling")
    ap.add_argument(
        "--max_cells_per_bucket",
        type=int,
        default=1000,
        help="Maximum cells per (cell_type, donor_id) bucket when metadata columns are present",
    )
    ap.add_argument("--n_programs", type=int, default=20, help="NMF factor count (number of programs)")
    ap.add_argument("--max_iter", type=int, default=500, help="NMF max iterations")
    ap.add_argument("--random_seed", type=int, default=1, help="Random seed")
    return ap.parse_args()


def main() -> int:
    args = _parse_args()
    np, pd, NMF = _load_optional_deps()

    matrix_path = Path(args.matrix_tsv)
    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    matrix_df = pd.read_csv(matrix_path, sep="\t")
    if args.cell_id_column not in matrix_df.columns:
        raise SystemExit(
            f"matrix_tsv is missing required cell ID column '{args.cell_id_column}'. "
            "Provide --cell_id_column to match your matrix."
        )
    cell_ids = matrix_df[args.cell_id_column].astype(str)
    expr_df = matrix_df.drop(columns=[args.cell_id_column])
    if expr_df.empty:
        raise SystemExit("matrix_tsv has no gene columns after removing cell_id column.")

    expr_numeric = expr_df.apply(pd.to_numeric, errors="coerce")
    n_non_numeric = int(expr_numeric.isna().sum().sum())
    if n_non_numeric > 0:
        print(
            f"warning: {n_non_numeric} non-numeric matrix values were coerced to 0.",
            file=sys.stderr,
        )
    expr_numeric = expr_numeric.fillna(0.0)

    sample_df = pd.DataFrame({args.cell_id_column: cell_ids})
    if args.metadata_tsv:
        meta_df = pd.read_csv(args.metadata_tsv, sep="\t")
        if args.cell_id_column not in meta_df.columns:
            raise SystemExit(
                f"metadata_tsv is missing required cell ID column '{args.cell_id_column}'. "
                "Provide --cell_id_column to match your metadata."
            )
        keep_cols = [args.cell_id_column]
        for col in (args.cell_type_column, args.donor_id_column):
            if col in meta_df.columns and col not in keep_cols:
                keep_cols.append(col)
        sample_df = sample_df.merge(meta_df[keep_cols], on=args.cell_id_column, how="left")

    bucket_cols: list[str] = []
    for col in (args.cell_type_column, args.donor_id_column):
        if col in sample_df.columns:
            bucket_cols.append(col)
    if bucket_cols:
        sample_df["_bucket"] = sample_df[bucket_cols].fillna("unknown").astype(str).agg("|".join, axis=1)
    else:
        sample_df["_bucket"] = "__all__"

    rng = np.random.default_rng(int(args.random_seed))
    selected_idx: list[int] = []
    for bucket, idx_arr in sample_df.groupby("_bucket").indices.items():
        idx_list = list(idx_arr)
        if len(idx_list) > int(args.max_cells_per_bucket):
            idx_list = list(rng.choice(idx_list, size=int(args.max_cells_per_bucket), replace=False))
        selected_idx.extend(idx_list)

    if len(selected_idx) > int(args.max_cells_total):
        selected_idx = list(rng.choice(selected_idx, size=int(args.max_cells_total), replace=False))

    selected_idx = sorted(set(int(i) for i in selected_idx))
    X = expr_numeric.iloc[selected_idx].to_numpy(dtype=float)
    gene_names = [str(c) for c in expr_numeric.columns]

    if X.size == 0:
        raise SystemExit("No cells remained after downsampling; increase --max_cells_total or --max_cells_per_bucket.")

    n_negative = int((X < 0).sum())
    if n_negative > 0:
        print(
            f"warning: found {n_negative} negative expression values; clipping to 0 for NMF.",
            file=sys.stderr,
        )
        X = np.clip(X, a_min=0.0, a_max=None)

    n_cells, n_genes = X.shape
    n_programs = int(args.n_programs)
    if n_programs <= 0:
        raise SystemExit("--n_programs must be > 0")
    if n_programs > min(n_cells, n_genes):
        raise SystemExit(
            f"--n_programs={n_programs} exceeds min(n_cells, n_genes)={min(n_cells, n_genes)} after downsampling."
        )

    model = NMF(
        n_components=n_programs,
        init="nndsvda",
        random_state=int(args.random_seed),
        max_iter=int(args.max_iter),
    )
    _ = model.fit_transform(X)
    H = model.components_  # shape: (n_programs, n_genes)

    out_df = pd.DataFrame({"gene_id": gene_names})
    for i in range(n_programs):
        out_df[f"program_{i + 1}"] = H[i, :]
    out_df.to_csv(out_path, sep="\t", index=False)

    print(
        "wrote loadings: "
        f"{out_path} (n_cells_used={n_cells}, n_genes={n_genes}, n_programs={n_programs})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
