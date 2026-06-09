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
    candidate = table_path.with_name(f"{table_path.stem}.provenance_graph.json")
    return str(candidate) if candidate.exists() else None


def _read_rows(args) -> list[dict[str, object]]:
    path = Path(args.table_tsv)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        required = [args.term_column, args.gene_symbol_column, args.score_column]
        missing = [name for name in required if name not in fieldnames]
        if missing:
            raise ValueError(f"Input TSV is missing required columns: {', '.join(missing)}")
        if args.gene_id_column and args.gene_id_column not in fieldnames:
            raise ValueError(f"Input TSV is missing gene_id_column: {args.gene_id_column}")
        if args.sign_column and args.sign_column not in fieldnames:
            raise ValueError(f"Input TSV is missing sign_column: {args.sign_column}")

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
            score = float(row.get(args.score_column, "0") or 0.0)
            if args.sign_column:
                sign_value = float(row.get(args.sign_column, "0") or 0.0)
            else:
                sign_value = 1.0 if score > 0 else (-1.0 if score < 0 else 0.0)
            direction = "up" if sign_value > 0 else ("dn" if sign_value < 0 else "")
            if not direction:
                continue
            out.append(
                {
                    "term": term,
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "score": abs(score),
                    "signed_score": score,
                    "sign": sign_value,
                    "direction": direction,
                }
            )
    return out


def _write_full_tables(out_dir: Path, rows: list[dict[str, object]]) -> None:
    sorted_rows = sorted(
        rows,
        key=lambda r: (str(r["term"]), str(r["direction"]), -float(r["score"]), str(r["gene_symbol"]), str(r["gene_id"])),
    )
    for index, row in enumerate(sorted_rows, start=1):
        row["rank"] = index
    fieldnames = ["term", "direction", "gene_id", "gene_symbol", "score", "signed_score", "sign", "rank"]
    for name in ["geneset.tsv", "geneset.full.tsv"]:
        path = out_dir / name
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(sorted_rows)


def _build_gene_sets_grouped_rows(args, rows: list[dict[str, object]]) -> tuple[list[tuple[str, list[str]]], list[dict[str, object]]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["term"]), str(row["direction"]))].append(row)

    gene_sets: list[tuple[str, list[str]]] = []
    summary_rows: list[dict[str, object]] = []
    label_map = {"up": "up", "dn": "dn"}
    if args.gmt_signed_labels == "pos_neg":
        label_map = {"up": "pos", "dn": "neg"}
    elif args.gmt_signed_labels == "Up_Down":
        label_map = {"up": "Up", "dn": "Down"}
    for (term, direction), group_rows in sorted(grouped.items()):
        genes = choose_gene_tokens(
            group_rows,
            prefer_symbol=bool(args.gmt_prefer_symbol),
            require_symbol=bool(args.gmt_require_symbol),
        )
        if len(genes) < int(args.gmt_min_genes) and not bool(args.emit_small_gene_sets):
            continue
        set_name = f"{term}{args.gmt_name_separator}{label_map[direction]}"
        gene_sets.append((set_name, genes))
        summary_rows.append(
            {
                "term": term,
                "direction": direction,
                "set_name": set_name,
                "gene_count": len(genes),
            }
        )
    return gene_sets, summary_rows


def _build_gene_sets_ternary_matrix_notebook(args, rows: list[dict[str, object]]) -> tuple[list[tuple[str, list[str]]], list[dict[str, object]]]:
    label_map = {"up": "up", "dn": "dn"}
    if args.gmt_signed_labels == "pos_neg":
        label_map = {"up": "pos", "dn": "neg"}
    elif args.gmt_signed_labels == "Up_Down":
        label_map = {"up": "Up", "dn": "Down"}

    term_gene_sign: dict[str, dict[str, float]] = defaultdict(dict)
    for row in rows:
        gene_symbol = str(row.get("gene_symbol", "")).strip()
        gene_id = str(row.get("gene_id", "")).strip()
        if bool(args.gmt_require_symbol):
            gene_token = gene_symbol
        elif bool(args.gmt_prefer_symbol):
            gene_token = gene_symbol or gene_id
        else:
            gene_token = gene_id or gene_symbol
        if not gene_token:
            continue
        term = str(row["term"])
        sign = float(row["sign"])
        existing = term_gene_sign[term].get(gene_token)
        if existing is None or sign > existing:
            term_gene_sign[term][gene_token] = sign

    gene_sets: list[tuple[str, list[str]]] = []
    summary_rows: list[dict[str, object]] = []
    for term in sorted(term_gene_sign):
        genes_by_term = term_gene_sign[term]
        up_genes = sorted(gene for gene, sign in genes_by_term.items() if sign > 0)
        dn_genes = sorted(gene for gene, sign in genes_by_term.items() if sign < 0)
        for direction, genes in [("up", up_genes), ("dn", dn_genes)]:
            if len(genes) < int(args.gmt_min_genes) and not bool(args.emit_small_gene_sets):
                continue
            set_name = f"{term}{args.gmt_name_separator}{label_map[direction]}"
            gene_sets.append((set_name, genes))
            summary_rows.append(
                {
                    "term": term,
                    "direction": direction,
                    "set_name": set_name,
                    "gene_count": len(genes),
                }
            )
    return gene_sets, summary_rows


def _build_gene_sets(args, rows: list[dict[str, object]]) -> tuple[list[tuple[str, list[str]]], list[dict[str, object]]]:
    emit_mode = str(getattr(args, "emit_mode", "grouped_rows") or "grouped_rows").strip()
    if emit_mode == "ternary_matrix_notebook":
        return _build_gene_sets_ternary_matrix_notebook(args, rows)
    return _build_gene_sets_grouped_rows(args, rows)


def run(args) -> dict[str, object]:
    activate_runtime_context("signed_term_gene", getattr(args, "provenance_overlay_json", None))
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
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["term", "direction", "set_name", "gene_count"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    meta = make_metadata(
        converter_name="signed_term_gene",
        parameters={
            "term_column": args.term_column,
            "gene_id_column": args.gene_id_column,
            "gene_symbol_column": args.gene_symbol_column,
            "score_column": args.score_column,
            "sign_column": args.sign_column,
            "emit_mode": args.emit_mode,
            "gmt_name_separator": args.gmt_name_separator,
            "gmt_signed_labels": args.gmt_signed_labels,
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
            "weight_type": "signed",
            "normalization": {"method": "none", "target_sum": None},
            "aggregation": "ternary_matrix_per_term" if args.emit_mode == "ternary_matrix_notebook" else "group_by_term_and_direction",
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
