from __future__ import annotations

import csv
from pathlib import Path

from omics2geneset.core.gmt import (
    build_gmt_sets_from_rows,
    parse_int_list_csv,
    parse_mass_list_csv,
    parse_str_list_csv,
    resolve_gmt_out_path,
    write_gmt,
)
from omics2geneset.core.atac_programs import (
    PROGRAM_DISTAL_ACTIVITY,
    PROGRAM_ENHANCER_BIAS,
    PROGRAM_LINKED_ACTIVITY,
    PROGRAM_PROMOTER_ACTIVITY,
    enhancer_bias_scores,
    mask_peak_weights,
    promoter_peak_indices,
    resolve_program_methods,
)
from omics2geneset.core.metadata import input_file_record, make_metadata, write_metadata
from omics2geneset.core.peak_to_gene import (
    link_distance_decay,
    link_external_regions,
    link_nearest_tss,
    link_promoter_overlap,
)
from omics2geneset.core.scoring import score_genes
from omics2geneset.core.selection import (
    global_l1_weights,
    ranked_gene_ids,
    select_quantile,
    select_threshold,
    select_top_k,
    within_set_l1_weights,
)
from omics2geneset.io.bed import read_bed, read_peak_weights_tsv
from omics2geneset.io.gtf import read_genes_from_gtf
from omics2geneset.io.region_gene_links import read_region_gene_links_tsv


_LINK_METHODS = ("promoter_overlap", "nearest_tss", "distance_decay")
_EXTERNAL_LINK_METHOD = "external"


def _arg(args, name: str, default):
    return getattr(args, name, default)


def _extract_peak_weights(peaks: list[dict[str, object]], peaks_weight_column: int | None, peak_weights_tsv: str | None) -> list[float]:
    if peak_weights_tsv:
        m = read_peak_weights_tsv(peak_weights_tsv)
        out: list[float] = []
        for p in peaks:
            key = (str(p["chrom"]), int(p["start"]), int(p["end"]))
            if key not in m:
                raise ValueError(f"Peak missing in peak_weights_tsv: {key}")
            out.append(float(m[key]))
        return out
    if peaks_weight_column is None:
        peaks_weight_column = 5
    idx = peaks_weight_column - 1
    out = []
    for p in peaks:
        cols = p.get("columns", [])
        if idx >= len(cols):
            raise ValueError("peaks_weight_column out of range for peaks file")
        out.append(float(cols[idx]))
    return out


def _resolve_link_methods(args) -> list[str]:
    link_expr = str(_arg(args, "link_method", "all")).strip()
    if not link_expr:
        link_expr = "all"
    tokens = [tok.strip() for tok in link_expr.split(",") if tok.strip()]
    if not tokens:
        tokens = ["all"]

    out: list[str] = []
    seen: set[str] = set()
    for tok in tokens:
        methods = list(_LINK_METHODS) if tok == "all" else [tok]
        for method in methods:
            if method not in _LINK_METHODS and method != _EXTERNAL_LINK_METHOD:
                raise ValueError(f"Unsupported link_method token: {method}")
            if method not in seen:
                out.append(method)
                seen.add(method)
    return out


def _resolve_max_distance(args, link_method: str) -> int:
    if link_method == _EXTERNAL_LINK_METHOD:
        return 0
    max_distance_bp = _arg(args, "max_distance_bp", None)
    if max_distance_bp is not None:
        return int(max_distance_bp)
    if link_method == "nearest_tss":
        return 100000
    if link_method == "distance_decay":
        return 500000
    return 100000


def _resolved_parameters(args) -> dict[str, object]:
    gmt_topk_list = parse_int_list_csv(str(_arg(args, "gmt_topk_list", "100,200,500")))
    gmt_mass_list = parse_mass_list_csv(str(_arg(args, "gmt_mass_list", "0.5,0.8,0.9")))
    gmt_biotype_allowlist = parse_str_list_csv(str(_arg(args, "gmt_biotype_allowlist", "protein_coding")))
    link_methods = _resolve_link_methods(args)
    program_methods = resolve_program_methods(
        "bulk",
        _arg(args, "program_preset", "default"),
        _arg(args, "program_methods", None),
    )
    max_distance = (
        _resolve_max_distance(args, link_methods[0])
        if len(link_methods) == 1
        else {method: _resolve_max_distance(args, method) for method in link_methods}
    )
    return {
        "link_method": _arg(args, "link_method", "all"),
        "link_methods_evaluated": link_methods,
        "primary_link_method": link_methods[0],
        "program_preset": _arg(args, "program_preset", "default"),
        "program_methods": _arg(args, "program_methods", None),
        "program_methods_evaluated": program_methods,
        "promoter_upstream_bp": _arg(args, "promoter_upstream_bp", 2000),
        "promoter_downstream_bp": _arg(args, "promoter_downstream_bp", 500),
        "max_distance_bp": max_distance,
        "decay_length_bp": _arg(args, "decay_length_bp", 50000),
        "max_genes_per_peak": _arg(args, "max_genes_per_peak", 5),
        "external_linking_enabled": _EXTERNAL_LINK_METHOD in link_methods,
        "peak_weight_transform": _arg(args, "peak_weight_transform", "abs"),
        "normalize": _arg(args, "normalize", "within_set_l1"),
        "selection_method": _arg(args, "select", "top_k"),
        "top_k": _arg(args, "top_k", 200),
        "quantile": _arg(args, "quantile", 0.01),
        "min_score": _arg(args, "min_score", 0.0),
        "emit_full": bool(_arg(args, "emit_full", True)),
        "emit_gmt": bool(_arg(args, "emit_gmt", True)),
        "gmt_prefer_symbol": bool(_arg(args, "gmt_prefer_symbol", True)),
        "gmt_require_symbol": bool(_arg(args, "gmt_require_symbol", True)),
        "gmt_biotype_allowlist": gmt_biotype_allowlist,
        "gmt_min_genes": int(_arg(args, "gmt_min_genes", 100)),
        "gmt_max_genes": int(_arg(args, "gmt_max_genes", 500)),
        "gmt_topk_list": gmt_topk_list,
        "gmt_mass_list": gmt_mass_list,
        "gmt_split_signed": bool(_arg(args, "gmt_split_signed", False)),
        "aggregation": "sum",
    }


def _link(peaks: list[dict[str, object]], genes, args, link_method: str):
    if link_method == _EXTERNAL_LINK_METHOD:
        region_gene_links = _arg(args, "_region_gene_links_resolved", None)
        if region_gene_links is None:
            raise ValueError("internal error: missing resolved external region-gene links")
        return link_external_regions(peaks, region_gene_links)
    max_distance_bp = _resolve_max_distance(args, link_method)
    if link_method == "promoter_overlap":
        return link_promoter_overlap(
            peaks,
            genes,
            _arg(args, "promoter_upstream_bp", 2000),
            _arg(args, "promoter_downstream_bp", 500),
        )
    if link_method == "nearest_tss":
        return link_nearest_tss(peaks, genes, max_distance_bp)
    if link_method == "distance_decay":
        return link_distance_decay(
            peaks,
            genes,
            max_distance_bp,
            _arg(args, "decay_length_bp", 50000),
            _arg(args, "max_genes_per_peak", 5),
        )
    raise ValueError(f"Unsupported link_method: {link_method}")


def _strip_version(gene_id: str) -> str:
    return gene_id.split(".", 1)[0]


def _resolve_external_link_gene_ids(
    region_gene_links: list[dict[str, object]],
    genes,
) -> tuple[list[dict[str, object]], int]:
    full_map = {str(g.gene_id): str(g.gene_id) for g in genes}
    nover_map: dict[str, str | None] = {}
    for g in genes:
        k = _strip_version(str(g.gene_id))
        if k in nover_map and nover_map[k] != str(g.gene_id):
            nover_map[k] = None
        else:
            nover_map[k] = str(g.gene_id)

    resolved: list[dict[str, object]] = []
    unresolved = 0
    for row in region_gene_links:
        gid = str(row["gene_id"])
        mapped = full_map.get(gid)
        if mapped is None:
            mapped = nover_map.get(_strip_version(gid))
        if not mapped:
            unresolved += 1
            continue
        resolved.append({**row, "gene_id": mapped})
    return resolved, unresolved


def _prepare_external_links(args, genes) -> tuple[list[dict[str, object]], int]:
    path = _arg(args, "region_gene_links_tsv", None)
    if not path:
        raise ValueError(
            "link_method includes 'external' but --region_gene_links_tsv was not provided"
        )
    raw = read_region_gene_links_tsv(path)
    resolved, unresolved = _resolve_external_link_gene_ids(raw, genes)
    if not resolved:
        raise ValueError(
            "No external region-gene links mapped to GTF gene_ids; check gene_id namespace/version"
        )
    return resolved, unresolved


def _select_gene_ids(scores: dict[str, float], args) -> list[str]:
    select = _arg(args, "select", "top_k")
    if select == "none":
        return ranked_gene_ids(scores)
    if select == "top_k":
        return select_top_k(scores, int(_arg(args, "top_k", 200)))
    if select == "quantile":
        return select_quantile(scores, float(_arg(args, "quantile", 0.01)))
    if select == "threshold":
        return select_threshold(scores, float(_arg(args, "min_score", 0.0)))
    raise ValueError(f"Unsupported selection method: {select}")


def _selected_weights(full_scores: dict[str, float], selected_gene_ids: list[str], normalize: str) -> dict[str, float]:
    if normalize == "none":
        return {g: float(full_scores[g]) for g in selected_gene_ids}
    if normalize == "l1":
        global_weights = global_l1_weights(full_scores)
        return {g: float(global_weights.get(g, 0.0)) for g in selected_gene_ids}
    if normalize == "within_set_l1":
        return within_set_l1_weights(full_scores, selected_gene_ids)
    raise ValueError(f"Unsupported normalization method: {normalize}")


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["gene_id", "score", "rank"]
    if any("weight" in row for row in rows):
        fieldnames.append("weight")
    if any("gene_symbol" in row and row["gene_symbol"] is not None for row in rows):
        fieldnames.append("gene_symbol")
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _rows_from_scores(
    scores: dict[str, float],
    gene_symbol_by_id: dict[str, str | None],
    gene_biotype_by_id: dict[str, str | None],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for rank, gene_id in enumerate(ranked_gene_ids(scores), start=1):
        rows.append(
            {
                "gene_id": gene_id,
                "score": float(scores[gene_id]),
                "rank": rank,
                "gene_symbol": gene_symbol_by_id.get(gene_id),
                "gene_biotype": gene_biotype_by_id.get(gene_id),
            }
        )
    return rows


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    peaks = read_bed(args.peaks)
    genes = read_genes_from_gtf(args.gtf)
    peak_weights = _extract_peak_weights(peaks, args.peaks_weight_column, args.peak_weights_tsv)
    program_methods = resolve_program_methods(
        "bulk",
        _arg(args, "program_preset", "default"),
        _arg(args, "program_methods", None),
    )
    link_methods = _resolve_link_methods(args)
    primary_link_method = link_methods[0]

    external_links_unresolved = 0
    if _EXTERNAL_LINK_METHOD in link_methods:
        resolved_external_links, external_links_unresolved = _prepare_external_links(args, genes)
        setattr(args, "_region_gene_links_resolved", resolved_external_links)
    links_by_method = {method: _link(peaks, genes, args, method) for method in link_methods}

    peak_weight_transform = _arg(args, "peak_weight_transform", "abs")
    normalize = _arg(args, "normalize", "within_set_l1")
    emit_full = bool(_arg(args, "emit_full", True))
    emit_gmt = bool(_arg(args, "emit_gmt", True))
    gmt_prefer_symbol = bool(_arg(args, "gmt_prefer_symbol", True))
    gmt_require_symbol = bool(_arg(args, "gmt_require_symbol", True))
    gmt_biotype_allowlist = parse_str_list_csv(str(_arg(args, "gmt_biotype_allowlist", "protein_coding")))
    gmt_min_genes = int(_arg(args, "gmt_min_genes", 100))
    gmt_max_genes = int(_arg(args, "gmt_max_genes", 500))
    gmt_topk_list = parse_int_list_csv(str(_arg(args, "gmt_topk_list", "100,200,500")))
    gmt_mass_list = parse_mass_list_csv(str(_arg(args, "gmt_mass_list", "0.5,0.8,0.9")))
    gmt_split_signed = bool(_arg(args, "gmt_split_signed", False))
    gmt_out = _arg(args, "gmt_out", None)

    gene_symbol_by_id = {g.gene_id: g.gene_symbol for g in genes}
    gene_biotype_by_id = {g.gene_id: g.gene_biotype for g in genes}

    full_scores_by_method: dict[str, dict[str, float]] = {}
    for method in link_methods:
        raw_scores = score_genes(peak_weights, links_by_method[method], peak_weight_transform)
        full_scores_by_method[method] = {g: float(s) for g, s in raw_scores.items() if float(s) != 0.0}
    full_scores = full_scores_by_method[primary_link_method]
    selected_gene_ids = _select_gene_ids(full_scores, args)
    selected_weights = _selected_weights(full_scores, selected_gene_ids, normalize)

    selected_rows: list[dict[str, object]] = []
    for rank, gene_id in enumerate(selected_gene_ids, start=1):
        selected_rows.append(
            {
                "gene_id": gene_id,
                "score": float(full_scores[gene_id]),
                "weight": float(selected_weights.get(gene_id, 0.0)),
                "rank": rank,
                "gene_symbol": gene_symbol_by_id.get(gene_id),
                "gene_biotype": gene_biotype_by_id.get(gene_id),
            }
        )
    _write_rows(out_dir / "geneset.tsv", selected_rows)

    full_rows_by_method: dict[str, list[dict[str, object]]] = {}
    for method in link_methods:
        full_rows_by_method[method] = _rows_from_scores(
            full_scores_by_method[method],
            gene_symbol_by_id,
            gene_biotype_by_id,
        )
    full_rows = full_rows_by_method[primary_link_method]
    if emit_full:
        _write_rows(out_dir / "geneset.full.tsv", full_rows)

    additional_program_rows: dict[str, list[dict[str, object]]] = {}
    need_promoter = (
        PROGRAM_PROMOTER_ACTIVITY in program_methods
        or PROGRAM_DISTAL_ACTIVITY in program_methods
        or PROGRAM_ENHANCER_BIAS in program_methods
    )
    need_distal = PROGRAM_DISTAL_ACTIVITY in program_methods or PROGRAM_ENHANCER_BIAS in program_methods
    promoter_links = (
        links_by_method["promoter_overlap"]
        if need_promoter and "promoter_overlap" in links_by_method
        else (_link(peaks, genes, args, "promoter_overlap") if need_promoter else [])
    )
    distance_links = (
        links_by_method["distance_decay"]
        if need_distal and "distance_decay" in links_by_method
        else (_link(peaks, genes, args, "distance_decay") if need_distal else [])
    )
    promoter_indices = promoter_peak_indices(promoter_links) if promoter_links else set()

    promoter_scores: dict[str, float] = {}
    distal_scores: dict[str, float] = {}
    if PROGRAM_PROMOTER_ACTIVITY in program_methods or PROGRAM_ENHANCER_BIAS in program_methods:
        promoter_peak_weights = mask_peak_weights(peak_weights, include_indices=promoter_indices)
        promoter_raw = score_genes(promoter_peak_weights, promoter_links, peak_weight_transform)
        promoter_scores = {g: float(s) for g, s in promoter_raw.items() if float(s) != 0.0}
        if PROGRAM_PROMOTER_ACTIVITY in program_methods:
            additional_program_rows[PROGRAM_PROMOTER_ACTIVITY] = _rows_from_scores(
                promoter_scores,
                gene_symbol_by_id,
                gene_biotype_by_id,
            )

    if PROGRAM_DISTAL_ACTIVITY in program_methods or PROGRAM_ENHANCER_BIAS in program_methods:
        distal_peak_weights = mask_peak_weights(peak_weights, exclude_indices=promoter_indices)
        distal_raw = score_genes(distal_peak_weights, distance_links, peak_weight_transform)
        distal_scores = {g: float(s) for g, s in distal_raw.items() if float(s) != 0.0}
        if PROGRAM_DISTAL_ACTIVITY in program_methods:
            additional_program_rows[PROGRAM_DISTAL_ACTIVITY] = _rows_from_scores(
                distal_scores,
                gene_symbol_by_id,
                gene_biotype_by_id,
            )

    if PROGRAM_ENHANCER_BIAS in program_methods:
        bias_scores = enhancer_bias_scores(promoter_scores, distal_scores)
        bias_nonzero = {g: float(s) for g, s in bias_scores.items() if float(s) != 0.0}
        additional_program_rows[PROGRAM_ENHANCER_BIAS] = _rows_from_scores(
            bias_nonzero,
            gene_symbol_by_id,
            gene_biotype_by_id,
        )

    gmt_path = resolve_gmt_out_path(out_dir, gmt_out)
    gmt_sets: list[tuple[str, list[str]]] = []
    gmt_plans: list[dict[str, object]] = []
    if emit_gmt:
        if PROGRAM_LINKED_ACTIVITY in program_methods:
            for method in link_methods:
                method_sets, method_plans = build_gmt_sets_from_rows(
                    rows=full_rows_by_method[method],
                    base_name=f"atac_bulk__program={PROGRAM_LINKED_ACTIVITY}__link_method={method}",
                    prefer_symbol=gmt_prefer_symbol,
                    min_genes=gmt_min_genes,
                    max_genes=gmt_max_genes,
                    topk_list=gmt_topk_list,
                    mass_list=gmt_mass_list,
                    split_signed=gmt_split_signed,
                    require_symbol=gmt_require_symbol,
                    allowed_biotypes={b.lower() for b in gmt_biotype_allowlist} if gmt_biotype_allowlist else None,
                )
                for plan in method_plans:
                    params = dict(plan.get("parameters", {}))
                    params["program_method"] = PROGRAM_LINKED_ACTIVITY
                    params["link_method"] = method
                    plan["parameters"] = params
                gmt_sets.extend(method_sets)
                gmt_plans.extend(method_plans)

        for program_method, rows in additional_program_rows.items():
            method_sets, method_plans = build_gmt_sets_from_rows(
                rows=rows,
                base_name=f"atac_bulk__program={program_method}",
                prefer_symbol=gmt_prefer_symbol,
                min_genes=gmt_min_genes,
                max_genes=gmt_max_genes,
                topk_list=gmt_topk_list,
                mass_list=gmt_mass_list,
                split_signed=gmt_split_signed,
                require_symbol=gmt_require_symbol,
                allowed_biotypes={b.lower() for b in gmt_biotype_allowlist} if gmt_biotype_allowlist else None,
            )
            for plan in method_plans:
                params = dict(plan.get("parameters", {}))
                params["program_method"] = program_method
                plan["parameters"] = params
            gmt_sets.extend(method_sets)
            gmt_plans.extend(method_plans)
        write_gmt(gmt_sets, gmt_path)

    assigned_peaks = len({int(link["peak_index"]) for link in links_by_method[primary_link_method]})
    files = [input_file_record(args.peaks, "peaks"), input_file_record(args.gtf, "gtf")]
    if args.peak_weights_tsv:
        files.append(input_file_record(args.peak_weights_tsv, "peak_weights_tsv"))
    if _EXTERNAL_LINK_METHOD in link_methods:
        files.append(input_file_record(args.region_gene_links_tsv, "region_gene_links_tsv"))

    output_files = [{"path": str(out_dir / "geneset.tsv"), "role": "selected_program"}]
    if emit_full:
        output_files.append({"path": str(out_dir / "geneset.full.tsv"), "role": "full_scores"})
    if emit_gmt:
        output_files.append({"path": str(gmt_path), "role": "gmt"})

    params = _resolved_parameters(args)
    meta = make_metadata(
        converter_name="atac_bulk",
        parameters=params,
        data_type="atac_seq",
        assay="bulk",
        organism=args.organism,
        genome_build=args.genome_build,
        files=files,
        gene_annotation={
            "mode": "gtf",
            "gtf_path": str(args.gtf),
            "source": args.gtf_source or "user",
            "gene_id_field": "gene_id",
        },
        weights={
            "weight_type": "signed" if peak_weight_transform == "signed" else "nonnegative",
            "normalization": {
                "method": normalize,
                "target_sum": 1.0 if normalize == "within_set_l1" else None,
            },
            "aggregation": "sum",
        },
        summary={
            "n_input_features": len(peaks),
            "n_genes": len(selected_rows),
            "n_features_assigned": assigned_peaks,
            "fraction_features_assigned": assigned_peaks / len(peaks) if peaks else 0.0,
            "n_external_links_unresolved_gene_id": external_links_unresolved,
            "n_program_methods": len(program_methods),
        },
        program_extraction={
            "selection_method": _arg(args, "select", "top_k"),
            "selection_params": {
                "k": _arg(args, "top_k", 200),
                "quantile": _arg(args, "quantile", 0.01),
                "min_score": _arg(args, "min_score", 0.0),
            },
            "normalize": normalize,
            "n_selected_genes": len(selected_rows),
            "score_definition": "sum over peaks of transformed peak_weight times link_weight",
            "program_methods": program_methods,
            "primary_program_method": PROGRAM_LINKED_ACTIVITY,
        },
        output_files=output_files,
        gmt={
            "written": emit_gmt,
            "path": (
                str(gmt_path.relative_to(out_dir)) if gmt_path.is_relative_to(out_dir) else str(gmt_path)
            )
            if emit_gmt
            else None,
            "prefer_symbol": gmt_prefer_symbol,
            "require_symbol": gmt_require_symbol,
            "biotype_allowlist": gmt_biotype_allowlist,
            "min_genes": gmt_min_genes,
            "max_genes": gmt_max_genes,
            "plans": gmt_plans,
        },
    )
    write_metadata(out_dir / "geneset.meta.json", meta)

    return {
        "n_peaks": len(peaks),
        "n_genes": len(selected_rows),
        "out_dir": str(out_dir),
    }
