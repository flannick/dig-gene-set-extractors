from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from geneset_extractors.core.gmt import choose_gene_tokens, write_gmt
from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context


def _resolve_upstream_provenance_graph_path(table_tsv: str | Path) -> str | None:
    table_path = Path(table_tsv)
    if not table_path.exists():
        return None
    candidates = [table_path.with_name(f"{table_path.stem}.provenance_graph.json")]
    if table_path.stem.endswith("_prefixed"):
        base_stem = table_path.stem[: -len("_prefixed")]
        candidates.append(table_path.with_name(f"{base_stem}.provenance_graph.json"))
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return None


def _read_rows(args) -> list[dict[str, object]]:
    path = Path(args.table_tsv)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        required = [args.term_column, args.gene_symbol_column]
        missing = [name for name in required if name not in fieldnames]
        if missing:
            raise ValueError(f"Input TSV is missing required columns: {', '.join(missing)}")
        if args.gene_id_column and args.gene_id_column not in fieldnames:
            raise ValueError(f"Input TSV is missing gene_id_column: {args.gene_id_column}")
        if args.score_column and args.score_column not in fieldnames:
            raise ValueError(f"Input TSV is missing score_column: {args.score_column}")

        out: list[dict[str, object]] = []
        for row in reader:
            term = str(row.get(args.term_column, "")).strip()
            gene_symbol = str(row.get(args.gene_symbol_column, "")).strip()
            gene_id = (
                str(row.get(args.gene_id_column, "")).strip()
                if args.gene_id_column
                else gene_symbol
            )
            if not gene_id:
                gene_id = gene_symbol
            if not term or not gene_symbol or not gene_id:
                continue
            score = float(row.get(args.score_column, "1") or 1.0) if args.score_column else 1.0
            out.append(
                {
                    "term": term,
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "score": score,
                }
            )
    return out


def _write_full_tables(out_dir: Path, rows: list[dict[str, object]]) -> None:
    sorted_rows = sorted(
        rows,
        key=lambda r: (str(r["term"]), -float(r["score"]), str(r["gene_symbol"]), str(r["gene_id"])),
    )
    for index, row in enumerate(sorted_rows, start=1):
        row["rank"] = index
        row["weight"] = float(row["score"])
    fieldnames = ["term", "gene_id", "gene_symbol", "score", "weight", "rank"]
    for name in ["geneset.tsv", "geneset.full.tsv"]:
        path = out_dir / name
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(sorted_rows)


def _build_gene_sets(args, rows: list[dict[str, object]]) -> tuple[list[tuple[str, list[str]]], list[dict[str, object]]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["term"])].append(row)

    gene_sets: list[tuple[str, list[str]]] = []
    summary_rows: list[dict[str, object]] = []
    for term in sorted(grouped):
        genes = choose_gene_tokens(
            grouped[term],
            prefer_symbol=bool(args.gmt_prefer_symbol),
            require_symbol=bool(args.gmt_require_symbol),
        )
        if len(genes) < int(args.gmt_min_genes) and not bool(args.emit_small_gene_sets):
            continue
        gene_sets.append((term, genes))
        summary_rows.append({"term": term, "set_name": term, "gene_count": len(genes)})
    return gene_sets, summary_rows


def run(args) -> dict[str, object]:
    activate_runtime_context("unsigned_term_gene", getattr(args, "provenance_overlay_json", None))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = _read_rows(args)
    _write_full_tables(out_dir, rows)
    gene_sets, summary_rows = _build_gene_sets(args, rows)
    upstream_graph_path = _resolve_upstream_provenance_graph_path(args.table_tsv)
    if bool(args.emit_gmt):
        write_gmt(gene_sets, out_dir / "genesets.gmt", gmt_format=getattr(args, "gmt_format", "classic"))

    summary_path = out_dir / "signature_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["term", "set_name", "gene_count"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    meta = make_metadata(
        converter_name="unsigned_term_gene",
        parameters={
            "term_column": args.term_column,
            "gene_id_column": args.gene_id_column,
            "gene_symbol_column": args.gene_symbol_column,
            "score_column": args.score_column,
            "gmt_min_genes": args.gmt_min_genes,
            "emit_small_gene_sets": args.emit_small_gene_sets,
        },
        data_type="transcriptomics",
        assay="bulk",
        organism=args.organism,
        genome_build=args.genome_build,
        files=[input_file_record(args.table_tsv, "table_tsv")],
        gene_annotation={"mode": "provided", "source": "input_table", "gene_id_field": args.gene_id_column or "gene_symbol"},
        weights={
            "weight_type": "nonnegative",
            "normalization": {"method": "none", "target_sum": None},
            "aggregation": "group_by_term",
        },
        summary={
            "n_input_features": len(rows),
            "n_genes": len({str(row["gene_id"]) for row in rows}),
            "n_features_assigned": len(rows),
            "fraction_features_assigned": 1.0 if rows else 0.0,
            "n_sets_emitted": len(gene_sets),
        },
        upstream_provenance_graph_path=upstream_graph_path,
        output_files=[
            {"path": "genesets.gmt", "role": "gmt_library"},
            {"path": "geneset.tsv", "role": "selected_program"},
            {"path": "geneset.full.tsv", "role": "full_scores"},
            {"path": "signature_summary.tsv", "role": "signature_summary"},
            {"path": "geneset.meta.json", "role": "metadata_json"},
        ],
    )
    write_metadata(out_dir / "geneset.meta.json", meta)
    return {"n_peaks": len(rows), "n_genes": len({str(row['gene_id']) for row in rows}), "out_dir": str(out_dir)}
