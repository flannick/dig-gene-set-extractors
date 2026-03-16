from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys

from geneset_extractors.core.gmt import (
    build_gmt_sets_from_rows,
    parse_int_list_csv,
    parse_mass_list_csv,
    parse_str_list_csv,
    resolve_gmt_out_path,
    write_gmt,
)
from geneset_extractors.core.metadata import make_metadata, write_metadata
from geneset_extractors.core.qc import write_run_summary_files
from geneset_extractors.core.selection import (
    global_l1_weights,
    ranked_gene_ids,
    select_quantile,
    select_threshold,
    select_top_k,
    within_set_l1_weights,
)
from geneset_extractors.extractors.rnaseq.deg_scoring import (
    DEGRow,
    annotate_records_with_gtf,
    apply_exclude_filters,
    compile_exclude_gene_regexes,
    maybe_promote_gene_id_to_symbol,
    parse_float_soft,
    resolve_deg_columns,
    sanitize_name_component,
    score_deg_rows,
)


@dataclass
class DEGWorkflowConfig:
    converter_name: str
    out_dir: Path
    organism: str
    genome_build: str
    signature_name: str
    deg_tsv_label: str
    comparison_label: str | None
    gene_id_column: str
    gene_symbol_column: str | None
    stat_column: str | None
    logfc_column: str | None
    padj_column: str | None
    pvalue_column: str | None
    score_column: str | None
    score_mode: str
    padj_max: float | None
    pvalue_max: float | None
    min_abs_logfc: float | None
    neglog10p_cap: float
    neglog10p_eps: float
    duplicate_gene_policy: str
    exclude_gene_regex: list[str] | None
    disable_default_excludes: bool
    gtf: str | None
    gtf_gene_id_field: str
    gtf_source: str | None
    select: str
    top_k: int
    quantile: float
    min_score: float
    normalize: str
    emit_full: bool
    emit_gmt: bool
    gmt_out: str | None
    gmt_prefer_symbol: bool
    gmt_require_symbol: bool
    gmt_biotype_allowlist: str
    gmt_min_genes: int
    gmt_max_genes: int
    gmt_topk_list: str
    gmt_mass_list: str
    gmt_split_signed: bool
    gmt_emit_abs: bool
    gmt_source: str
    emit_small_gene_sets: bool
    warn_biotype_missing: bool


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["gene_id", "score"]
    if any("weight" in row for row in rows):
        fieldnames.append("weight")
    if any(str(row.get("gene_symbol", "")).strip() for row in rows):
        fieldnames.append("gene_symbol")
    if any(str(row.get("gene_biotype", "")).strip() for row in rows):
        fieldnames.append("gene_biotype")
    if any("rank" in row for row in rows):
        fieldnames.append("rank")
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _select_gene_ids(scores: dict[str, float], cfg: DEGWorkflowConfig) -> list[str]:
    if cfg.select == "none":
        return ranked_gene_ids(scores)
    if cfg.select == "top_k":
        return select_top_k(scores, int(cfg.top_k))
    if cfg.select == "quantile":
        return select_quantile(scores, float(cfg.quantile))
    if cfg.select == "threshold":
        return select_threshold(scores, float(cfg.min_score))
    raise ValueError(f"Unsupported selection method: {cfg.select}")


def _selected_weights(magnitude_scores: dict[str, float], selected_gene_ids: list[str], normalize: str) -> dict[str, float]:
    if normalize == "none":
        return {g: float(magnitude_scores.get(g, 0.0)) for g in selected_gene_ids}
    if normalize == "l1":
        global_weights = global_l1_weights(magnitude_scores)
        return {g: float(global_weights.get(g, 0.0)) for g in selected_gene_ids}
    if normalize == "within_set_l1":
        return within_set_l1_weights(magnitude_scores, selected_gene_ids)
    raise ValueError(f"Unsupported normalization method: {normalize}")


def _collect_skipped_gmt_outputs(gmt_diagnostics: list[dict[str, object]]) -> list[dict[str, object]]:
    skipped: list[dict[str, object]] = []
    for diag in gmt_diagnostics:
        code = str(diag.get("code", ""))
        if code not in {"small_gene_set_skipped", "no_positive_genes"}:
            continue
        skipped.append(
            {
                "reason": str(diag.get("reason", "")),
                "code": code,
                "n_genes": int(diag.get("n_genes", 0) or 0),
                "variant": str(diag.get("variant", "primary")),
            }
        )
    return skipped


def _warn_gmt_diagnostics(gmt_diagnostics: list[dict[str, object]]) -> None:
    for diag in gmt_diagnostics:
        code = str(diag.get("code", ""))
        base_name = str(diag.get("base_name", "unknown"))
        n_genes = int(diag.get("n_genes", 0) or 0)
        reason = str(diag.get("reason", "")).strip()
        suggestion = str(diag.get("suggestion", "")).strip()
        if code in {"small_gene_set_skipped", "no_positive_genes"}:
            print(f"warning: skipped GMT output for {base_name}; n_genes={n_genes}. {reason}", file=sys.stderr)
        elif code == "small_gene_set_emitted":
            print(f"warning: emitted small GMT output for {base_name}; n_genes={n_genes}. {reason}", file=sys.stderr)
        elif code == "marginal_gene_count":
            print(f"warning: marginal GMT signal for {base_name}; n_genes={n_genes}. {reason}", file=sys.stderr)
        elif code == "require_symbol_heavy_drop":
            print(
                f"warning: GMT symbol requirement dropped many rows for {base_name}. {reason}",
                file=sys.stderr,
            )
        else:
            continue
        if suggestion:
            print(f"warning:   {suggestion}", file=sys.stderr)


def _warn_symbol_availability(full_rows: list[dict[str, object]], cfg: DEGWorkflowConfig) -> None:
    if not cfg.gmt_prefer_symbol:
        return
    if not full_rows:
        return
    n_with_symbol = sum(1 for row in full_rows if str(row.get("gene_symbol", "")).strip())
    n_total = len(full_rows)
    if n_with_symbol == 0:
        if cfg.gmt_require_symbol:
            print(
                "warning: no gene_symbol values found and gmt_require_symbol=true; "
                "GMT tokenization may drop most or all genes. Provide --gtf to map symbols, "
                "or set --gmt_require_symbol false.",
                file=sys.stderr,
            )
        else:
            print(
                "warning: no gene_symbol values found; GMT token selection will fall back to gene_id.",
                file=sys.stderr,
            )
        return
    if cfg.gmt_require_symbol and n_with_symbol < max(1, int(0.2 * n_total)):
        print(
            "warning: gmt_require_symbol=true but most rows lack gene_symbol; many genes may be dropped from GMT output.",
            file=sys.stderr,
        )


def _warn_split_signed_no_negatives(rows: list[dict[str, object]], cfg: DEGWorkflowConfig) -> None:
    if not cfg.gmt_split_signed:
        return
    if not rows:
        return
    has_negative = any(float(row.get("score", 0.0)) < 0.0 for row in rows)
    if has_negative:
        return
    print(
        "warning: No negative scores detected; writing unsplit (pos-only) GMT set(s).",
        file=sys.stderr,
    )


def _warn_biotype_availability(rows: list[dict[str, object]], cfg: DEGWorkflowConfig, allowlist: list[str]) -> bool:
    if not allowlist:
        return False
    if not rows:
        return False
    missing = 0
    for row in rows:
        value = row.get("gene_biotype")
        if value is None or not str(value).strip():
            missing += 1
    fraction = float(missing) / float(len(rows))
    if fraction > 0.8:
        print(
            "warning: gmt_biotype_allowlist is set but gene_biotype values are missing for "
            f"{missing}/{len(rows)} rows; filter may be weakly informative.",
            file=sys.stderr,
        )
        return True
    return False


def run_deg_workflow(
    *,
    cfg: DEGWorkflowConfig,
    fieldnames: list[str],
    rows: list[DEGRow],
    input_files: list[dict[str, str]],
) -> dict[str, object]:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if cfg.gmt_source not in {"full", "selected"}:
        raise ValueError(f"Unsupported gmt_source: {cfg.gmt_source}")

    resolved_columns = resolve_deg_columns(
        fieldnames,
        gene_id_column=cfg.gene_id_column,
        gene_symbol_column=cfg.gene_symbol_column,
        stat_column=cfg.stat_column,
        logfc_column=cfg.logfc_column,
        padj_column=cfg.padj_column,
        pvalue_column=cfg.pvalue_column,
        score_column=cfg.score_column,
    )
    if cfg.score_mode == "custom_column":
        custom_col = str(resolved_columns.get("score_column") or "")
        if custom_col not in fieldnames:
            raise ValueError(
                "score_mode=custom_column requires --score_column with a valid table column. "
                f"Provided: {cfg.score_column!r}; available: {', '.join(fieldnames)}"
            )

    filtered_input_rows = rows
    n_rows_filtered_by_thresholds = 0
    if cfg.padj_max is not None:
        padj_column = resolved_columns["padj_column"]
        if not padj_column:
            raise ValueError("--padj_max requires a resolvable adjusted p-value column")
        before = len(filtered_input_rows)
        filtered_input_rows = [
            row
            for row in filtered_input_rows
            if (parse_float_soft(row.values.get(str(padj_column))) is not None)
            and float(parse_float_soft(row.values.get(str(padj_column))) or 0.0) <= float(cfg.padj_max)
        ]
        n_rows_filtered_by_thresholds += before - len(filtered_input_rows)
    if cfg.pvalue_max is not None:
        pvalue_column = resolved_columns["pvalue_column"]
        if not pvalue_column:
            raise ValueError("--pvalue_max requires a resolvable p-value column")
        before = len(filtered_input_rows)
        filtered_input_rows = [
            row
            for row in filtered_input_rows
            if (parse_float_soft(row.values.get(str(pvalue_column))) is not None)
            and float(parse_float_soft(row.values.get(str(pvalue_column))) or 0.0) <= float(cfg.pvalue_max)
        ]
        n_rows_filtered_by_thresholds += before - len(filtered_input_rows)
    if cfg.min_abs_logfc is not None:
        logfc_column = resolved_columns["logfc_column"]
        if not logfc_column:
            raise ValueError("--min_abs_logfc requires a resolvable logFC column")
        before = len(filtered_input_rows)
        filtered_input_rows = [
            row
            for row in filtered_input_rows
            if (parse_float_soft(row.values.get(str(logfc_column))) is not None)
            and abs(float(parse_float_soft(row.values.get(str(logfc_column))) or 0.0)) >= float(cfg.min_abs_logfc)
        ]
        n_rows_filtered_by_thresholds += before - len(filtered_input_rows)
    if not filtered_input_rows:
        raise ValueError("No DE rows remain after applying padj/pvalue/logFC row filters.")

    resolved_score_mode, records, skipped_rows, duplicate_info = score_deg_rows(
        filtered_input_rows,
        gene_id_column=str(resolved_columns["gene_id_column"]),
        gene_symbol_column=resolved_columns["gene_symbol_column"],
        stat_column=resolved_columns["stat_column"],
        logfc_column=resolved_columns["logfc_column"],
        padj_column=resolved_columns["padj_column"],
        pvalue_column=resolved_columns["pvalue_column"],
        score_column=resolved_columns["score_column"],
        score_mode=cfg.score_mode,
        neglog10p_eps=float(cfg.neglog10p_eps),
        neglog10p_cap=float(cfg.neglog10p_cap),
        duplicate_gene_policy=cfg.duplicate_gene_policy,
    )
    n_input_features = len(records)
    if int(duplicate_info.get("n_gene_ids_with_duplicates", 0)) > 0:
        examples = duplicate_info.get("top_examples", [])
        preview = ", ".join(
            f"{str(item.get('gene_id', ''))}(n={int(item.get('n_rows', 0))})"
            for item in examples[:3]
        )
        print(
            "warning: detected duplicate gene_id rows in DE input; "
            f"n_gene_ids_with_duplicates={int(duplicate_info.get('n_gene_ids_with_duplicates', 0))}, "
            f"policy={cfg.duplicate_gene_policy}. "
            + (f"Examples: {preview}" if preview else ""),
            file=sys.stderr,
        )

    if cfg.gtf:
        records = annotate_records_with_gtf(records, cfg.gtf, cfg.gtf_gene_id_field)

    symbol_promotion = maybe_promote_gene_id_to_symbol(records)
    records = symbol_promotion["records"]
    if bool(symbol_promotion.get("promoted", False)):
        print(
            "warning: promoted gene_id values to gene_symbol for rows missing symbols "
            f"(n_promoted={int(symbol_promotion.get('n_promoted', 0))}).",
            file=sys.stderr,
        )
    elif str(symbol_promotion.get("reason", "")) == "gene_id_looks_ensembl_like":
        print(
            "warning: gene_symbol values are mostly missing and gene_id appears Ensembl-like; "
            "symbols were not auto-promoted. Provide --gtf for symbol mapping or set "
            "--gmt_require_symbol false to allow gene_id tokens.",
            file=sys.stderr,
        )

    exclude_patterns = compile_exclude_gene_regexes(
        cfg.exclude_gene_regex,
        cfg.disable_default_excludes,
    )
    records, n_filtered = apply_exclude_filters(records, exclude_patterns)
    if not records:
        raise ValueError("No genes remain after applying score parsing and filters.")

    score_by_gene = {gid: float(rec["score"]) for gid, rec in records.items()}
    magnitude_by_gene = {gid: abs(float(score)) for gid, score in score_by_gene.items()}

    selected_gene_ids = _select_gene_ids(magnitude_by_gene, cfg)
    selected_gene_ids = sorted(
        selected_gene_ids,
        key=lambda gid: (-float(magnitude_by_gene.get(gid, 0.0)), str(gid)),
    )
    weights = _selected_weights(magnitude_by_gene, selected_gene_ids, cfg.normalize)

    selected_rows: list[dict[str, object]] = []
    for rank, gid in enumerate(selected_gene_ids, start=1):
        rec = records[gid]
        row: dict[str, object] = {
            "gene_id": gid,
            "score": float(score_by_gene[gid]),
            "weight": float(weights.get(gid, 0.0)),
            "rank": rank,
        }
        symbol_obj = rec.get("gene_symbol")
        symbol = "" if symbol_obj is None else str(symbol_obj).strip()
        if symbol:
            row["gene_symbol"] = symbol
        biotype_obj = rec.get("gene_biotype")
        biotype = "" if biotype_obj is None else str(biotype_obj).strip()
        if biotype:
            row["gene_biotype"] = biotype
        selected_rows.append(row)

    full_gene_ids = sorted(
        [gid for gid, score in score_by_gene.items() if float(score) != 0.0],
        key=lambda gid: (-float(abs(score_by_gene[gid])), str(gid)),
    )
    full_rows: list[dict[str, object]] = []
    for rank, gid in enumerate(full_gene_ids, start=1):
        rec = records[gid]
        row = {
            "gene_id": gid,
            "score": float(score_by_gene[gid]),
            "rank": rank,
        }
        symbol_obj = rec.get("gene_symbol")
        symbol = "" if symbol_obj is None else str(symbol_obj).strip()
        if symbol:
            row["gene_symbol"] = symbol
        biotype_obj = rec.get("gene_biotype")
        biotype = "" if biotype_obj is None else str(biotype_obj).strip()
        if biotype:
            row["gene_biotype"] = biotype
        full_rows.append(row)

    _write_rows(out_dir / "geneset.tsv", selected_rows)
    if cfg.emit_full:
        _write_rows(out_dir / "geneset.full.tsv", full_rows)

    gmt_sets: list[tuple[str, list[str]]] = []
    gmt_plans: list[dict[str, object]] = []
    gmt_diagnostics: list[dict[str, object]] = []
    gmt_path = resolve_gmt_out_path(out_dir, cfg.gmt_out)
    gmt_topk_list = parse_int_list_csv(str(cfg.gmt_topk_list))
    gmt_mass_list = parse_mass_list_csv(str(cfg.gmt_mass_list))
    gmt_biotype_allowlist = parse_str_list_csv(str(cfg.gmt_biotype_allowlist))
    gmt_rows = full_rows if cfg.gmt_source == "full" else selected_rows
    if cfg.emit_gmt:
        _warn_symbol_availability(gmt_rows, cfg)
        biotype_warning_emitted = (
            _warn_biotype_availability(gmt_rows, cfg, gmt_biotype_allowlist)
            if cfg.warn_biotype_missing
            else False
        )
        _warn_split_signed_no_negatives(gmt_rows, cfg)
        safe_signature = sanitize_name_component(cfg.signature_name)
        safe_score_mode = sanitize_name_component(resolved_score_mode)
        base_name = f"{cfg.converter_name}__signature={safe_signature}__score_mode={safe_score_mode}"
        if cfg.comparison_label:
            safe_comparison = sanitize_name_component(cfg.comparison_label)
            base_name = (
                f"{cfg.converter_name}__comparison={safe_comparison}__signature={safe_signature}"
                f"__score_mode={safe_score_mode}"
            )
        gmt_sets, gmt_plans = build_gmt_sets_from_rows(
            rows=gmt_rows,
            base_name=base_name,
            prefer_symbol=bool(cfg.gmt_prefer_symbol),
            min_genes=int(cfg.gmt_min_genes),
            max_genes=int(cfg.gmt_max_genes),
            topk_list=gmt_topk_list,
            mass_list=gmt_mass_list,
            split_signed=bool(cfg.gmt_split_signed),
            require_symbol=bool(cfg.gmt_require_symbol),
            allowed_biotypes={x.lower() for x in gmt_biotype_allowlist} if gmt_biotype_allowlist else None,
            emit_small_gene_sets=bool(cfg.emit_small_gene_sets),
            diagnostics=gmt_diagnostics,
            context={
                "converter": cfg.converter_name,
                "comparison": cfg.comparison_label or "",
            },
        )
        if cfg.gmt_emit_abs:
            abs_rows = [
                {
                    "gene_id": str(row.get("gene_id", "")),
                    "gene_symbol": row.get("gene_symbol"),
                    "gene_biotype": row.get("gene_biotype"),
                    "score": abs(float(row.get("score", 0.0))),
                    "rank": row.get("rank"),
                }
                for row in gmt_rows
            ]
            abs_sets, abs_plans = build_gmt_sets_from_rows(
                rows=abs_rows,
                base_name=f"{base_name}__abs",
                prefer_symbol=bool(cfg.gmt_prefer_symbol),
                min_genes=int(cfg.gmt_min_genes),
                max_genes=int(cfg.gmt_max_genes),
                topk_list=gmt_topk_list,
                mass_list=gmt_mass_list,
                split_signed=False,
                require_symbol=bool(cfg.gmt_require_symbol),
                allowed_biotypes={x.lower() for x in gmt_biotype_allowlist} if gmt_biotype_allowlist else None,
                emit_small_gene_sets=bool(cfg.emit_small_gene_sets),
                diagnostics=gmt_diagnostics,
                context={
                    "converter": cfg.converter_name,
                    "comparison": cfg.comparison_label or "",
                    "direction": "abs",
                },
            )
            gmt_sets.extend(abs_sets)
            gmt_plans.extend(abs_plans)
        write_gmt(gmt_sets, gmt_path)
    _warn_gmt_diagnostics(gmt_diagnostics)

    output_files = [{"path": str(out_dir / "geneset.tsv"), "role": "selected_program"}]
    if cfg.emit_full:
        output_files.append({"path": str(out_dir / "geneset.full.tsv"), "role": "full_scores"})
    if cfg.emit_gmt:
        output_files.append({"path": str(gmt_path), "role": "gmt"})

    run_summary_payload: dict[str, object] = {
        "converter": cfg.converter_name,
        "signature_name": cfg.signature_name,
        "comparison": cfg.comparison_label,
        "score_mode": resolved_score_mode,
        "n_input_features": n_input_features,
        "n_genes_after_filter": len(records),
        "n_genes_selected": len(selected_rows),
        "n_rows_skipped_unparseable": skipped_rows,
        "n_rows_filtered_by_thresholds": n_rows_filtered_by_thresholds,
        "n_genes_filtered_by_symbol_regex": n_filtered,
        "duplicate_gene_policy": cfg.duplicate_gene_policy,
        "n_gene_ids_with_duplicates": int(duplicate_info.get("n_gene_ids_with_duplicates", 0)),
        "n_duplicate_rows": int(duplicate_info.get("n_duplicate_rows", 0)),
        "duplicate_gene_examples": duplicate_info.get("top_examples", []),
        "symbol_promotion": {
            "promoted": bool(symbol_promotion.get("promoted", False)),
            "reason": str(symbol_promotion.get("reason", "")),
            "n_promoted": int(symbol_promotion.get("n_promoted", 0) or 0),
            "ensembl_like_fraction": float(symbol_promotion.get("ensembl_like_fraction", 0.0) or 0.0),
        },
        "requested_gmt_outputs": [{"name": p.get("name", "")} for p in gmt_plans],
        "skipped_gmt_outputs": _collect_skipped_gmt_outputs(gmt_diagnostics),
    }
    run_summary_json_path, run_summary_txt_path = write_run_summary_files(out_dir, run_summary_payload)
    output_files.append({"path": str(run_summary_json_path), "role": "run_summary_json"})
    output_files.append({"path": str(run_summary_txt_path), "role": "run_summary_text"})

    params: dict[str, object] = {
        "signature_name": cfg.signature_name,
        "comparison_label": cfg.comparison_label,
        "columns": {
            "gene_id_column": resolved_columns["gene_id_column"],
            "gene_symbol_column": resolved_columns["gene_symbol_column"],
            "stat_column": resolved_columns["stat_column"],
            "logfc_column": resolved_columns["logfc_column"],
            "padj_column": resolved_columns["padj_column"],
            "pvalue_column": resolved_columns["pvalue_column"],
            "score_column": resolved_columns["score_column"],
        },
        "score_mode": resolved_score_mode,
        "score_mode_requested": cfg.score_mode,
        "padj_max": cfg.padj_max,
        "pvalue_max": cfg.pvalue_max,
        "min_abs_logfc": cfg.min_abs_logfc,
        "duplicate_gene_policy": cfg.duplicate_gene_policy,
        "neglog10p_cap": cfg.neglog10p_cap,
        "neglog10p_eps": cfg.neglog10p_eps,
        "selection_method": cfg.select,
        "top_k": cfg.top_k,
        "quantile": cfg.quantile,
        "min_score": cfg.min_score,
        "normalize": cfg.normalize,
        "exclude_gene_regex": cfg.exclude_gene_regex or [],
        "disable_default_excludes": cfg.disable_default_excludes,
        "emit_full": cfg.emit_full,
        "emit_gmt": cfg.emit_gmt,
        "gmt_split_signed": cfg.gmt_split_signed,
        "gmt_emit_abs": cfg.gmt_emit_abs,
        "gmt_topk_list": gmt_topk_list,
        "gmt_mass_list": gmt_mass_list,
        "gmt_source": cfg.gmt_source,
        "gmt_min_genes": cfg.gmt_min_genes,
        "gmt_max_genes": cfg.gmt_max_genes,
        "gmt_prefer_symbol": cfg.gmt_prefer_symbol,
        "gmt_require_symbol": cfg.gmt_require_symbol,
        "gmt_biotype_allowlist": gmt_biotype_allowlist,
        "emit_small_gene_sets": cfg.emit_small_gene_sets,
    }

    meta = make_metadata(
        converter_name=cfg.converter_name,
        parameters=params,
        data_type="rna_seq",
        assay="bulk",
        organism=cfg.organism,
        genome_build=cfg.genome_build,
        files=input_files,
        gene_annotation=(
            {
                "mode": "gtf",
                "gtf_path": str(cfg.gtf),
                "source": cfg.gtf_source or "user",
                "gene_id_field": cfg.gtf_gene_id_field,
            }
            if cfg.gtf
            else {
                "mode": "none",
                "source": "none",
                "gene_id_field": str(resolved_columns["gene_id_column"]),
            }
        ),
        weights={
            "weight_type": "signed",
            "normalization": {
                "method": cfg.normalize,
                "target_sum": 1.0 if cfg.normalize in {"within_set_l1", "l1"} else None,
            },
            "aggregation": cfg.duplicate_gene_policy,
        },
        summary={
            "n_input_features": int(n_input_features),
            "n_genes": int(len(selected_rows)),
            "n_features_assigned": int(n_input_features),
            "fraction_features_assigned": 1.0 if n_input_features else 0.0,
            "n_rows_skipped_unparseable": int(skipped_rows),
            "n_rows_filtered_by_thresholds": int(n_rows_filtered_by_thresholds),
            "n_genes_filtered_by_symbol_regex": int(n_filtered),
            "n_gene_ids_with_duplicates": int(duplicate_info.get("n_gene_ids_with_duplicates", 0)),
            "n_duplicate_rows": int(duplicate_info.get("n_duplicate_rows", 0)),
            "symbol_promotion": {
                "promoted": bool(symbol_promotion.get("promoted", False)),
                "reason": str(symbol_promotion.get("reason", "")),
                "n_promoted": int(symbol_promotion.get("n_promoted", 0) or 0),
                "ensembl_like_fraction": float(symbol_promotion.get("ensembl_like_fraction", 0.0) or 0.0),
            },
        },
        program_extraction={
            "selection_method": cfg.select,
            "selection_params": {
                "k": cfg.top_k,
                "quantile": cfg.quantile,
                "min_score": cfg.min_score,
            },
            "normalize": cfg.normalize,
            "n_selected_genes": len(selected_rows),
            "score_definition": "signed differential expression score aggregated by gene_id; selection/weights on abs(score)",
        },
        output_files=output_files,
        gmt={
            "written": bool(cfg.emit_gmt),
            "path": (
                str(gmt_path.relative_to(out_dir)) if cfg.emit_gmt and gmt_path.is_relative_to(out_dir) else str(gmt_path)
            )
            if cfg.emit_gmt
            else None,
            "prefer_symbol": bool(cfg.gmt_prefer_symbol),
            "require_symbol": bool(cfg.gmt_require_symbol),
            "biotype_allowlist": gmt_biotype_allowlist,
            "source": cfg.gmt_source,
            "emit_abs": bool(cfg.gmt_emit_abs),
            "min_genes": int(cfg.gmt_min_genes),
            "max_genes": int(cfg.gmt_max_genes),
            "emit_small_gene_sets": bool(cfg.emit_small_gene_sets),
            "requested_outputs": [{"name": p.get("name", "")} for p in gmt_plans],
            "emitted_outputs": [{"name": p.get("name", "")} for p in gmt_plans],
            "skipped_outputs": _collect_skipped_gmt_outputs(gmt_diagnostics),
            "diagnostics": gmt_diagnostics,
            "plans": gmt_plans,
        },
    )
    write_metadata(out_dir / "geneset.meta.json", meta)

    return {
        "resolved_score_mode": resolved_score_mode,
        "n_input_features": n_input_features,
        "n_genes_selected": len(selected_rows),
        "n_rows_skipped_unparseable": skipped_rows,
        "n_rows_filtered_by_thresholds": n_rows_filtered_by_thresholds,
        "n_genes_filtered_by_symbol_regex": n_filtered,
        "biotype_warning_emitted": bool(biotype_warning_emitted if cfg.emit_gmt else False),
        "gmt_sets": gmt_sets,
        "gmt_plans": gmt_plans,
    }
