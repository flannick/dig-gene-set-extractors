from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import re

from omics2geneset.core.gmt import (
    build_gmt_sets_from_rows,
    parse_int_list_csv,
    parse_mass_list_csv,
    parse_str_list_csv,
    resolve_gmt_out_path,
    write_gmt,
)
from omics2geneset.core.metadata import input_file_record, make_metadata, write_metadata
from omics2geneset.core.peak_to_gene import link_distance_decay, link_nearest_tss, link_promoter_overlap
from omics2geneset.core.scoring import score_genes
from omics2geneset.core.selection import (
    global_l1_weights,
    ranked_gene_ids,
    select_quantile,
    select_threshold,
    select_top_k,
    within_set_l1_weights,
)
from omics2geneset.io.gtf import read_genes_from_gtf
from omics2geneset.io.mtx_10x import (
    make_group_indices,
    read_10x_matrix_dir,
    read_groups_tsv,
    summarize_group_vs_rest,
    summarize_peaks,
    summarize_peaks_by_group,
)


_LINK_METHODS = ("promoter_overlap", "nearest_tss", "distance_decay")


def _arg(args, name: str, default):
    return getattr(args, name, default)


def _resolve_link_methods(args) -> list[str]:
    link_method = str(_arg(args, "link_method", "all"))
    if link_method == "all":
        return list(_LINK_METHODS)
    if link_method not in _LINK_METHODS:
        raise ValueError(f"Unsupported link_method: {link_method}")
    return [link_method]


def _resolve_max_distance(args, link_method: str) -> int:
    max_distance_bp = _arg(args, "max_distance_bp", None)
    if max_distance_bp is not None:
        return int(max_distance_bp)
    if link_method == "nearest_tss":
        return 100000
    if link_method == "distance_decay":
        return 500000
    return 100000


def _resolve_contrast(args) -> str:
    contrast = _arg(args, "contrast", None)
    if contrast is not None:
        return str(contrast)
    if _arg(args, "groups_tsv", None):
        return "group_vs_rest"
    return "none"


def _resolve_pseudocount(args) -> float:
    contrast_pseudocount = _arg(args, "contrast_pseudocount", None)
    if contrast_pseudocount is not None:
        return float(contrast_pseudocount)
    if _arg(args, "peak_summary", "sum_counts") == "frac_cells_nonzero":
        return 1e-3
    return 1.0


def _resolved_parameters(args, group_name: str | None = None) -> dict[str, object]:
    contrast = _resolve_contrast(args)
    gmt_topk_list = parse_int_list_csv(str(_arg(args, "gmt_topk_list", "100,200,500")))
    gmt_mass_list = parse_mass_list_csv(str(_arg(args, "gmt_mass_list", "0.5,0.8,0.9")))
    gmt_biotype_allowlist = parse_str_list_csv(str(_arg(args, "gmt_biotype_allowlist", "protein_coding")))
    link_methods = _resolve_link_methods(args)
    max_distance = (
        _resolve_max_distance(args, link_methods[0])
        if len(link_methods) == 1
        else {method: _resolve_max_distance(args, method) for method in link_methods}
    )
    params: dict[str, object] = {
        "link_method": _arg(args, "link_method", "all"),
        "link_methods_evaluated": link_methods,
        "primary_link_method": link_methods[0],
        "promoter_upstream_bp": _arg(args, "promoter_upstream_bp", 2000),
        "promoter_downstream_bp": _arg(args, "promoter_downstream_bp", 500),
        "max_distance_bp": max_distance,
        "decay_length_bp": _arg(args, "decay_length_bp", 50000),
        "max_genes_per_peak": _arg(args, "max_genes_per_peak", 5),
        "peak_summary": _arg(args, "peak_summary", "sum_counts"),
        "peak_weight_transform": _arg(args, "peak_weight_transform", "positive"),
        "contrast": contrast,
        "contrast_metric": _arg(args, "contrast_metric", "log2fc"),
        "contrast_pseudocount": _resolve_pseudocount(args),
        "selection_method": _arg(args, "select", "top_k"),
        "top_k": _arg(args, "top_k", 200),
        "quantile": _arg(args, "quantile", 0.01),
        "min_score": _arg(args, "min_score", 0.0),
        "emit_full": bool(_arg(args, "emit_full", True)),
        "normalize": _arg(args, "normalize", "within_set_l1"),
        "emit_gmt": bool(_arg(args, "emit_gmt", True)),
        "gmt_prefer_symbol": bool(_arg(args, "gmt_prefer_symbol", True)),
        "gmt_require_symbol": bool(_arg(args, "gmt_require_symbol", True)),
        "gmt_biotype_allowlist": gmt_biotype_allowlist,
        "gmt_min_genes": int(_arg(args, "gmt_min_genes", 100)),
        "gmt_max_genes": int(_arg(args, "gmt_max_genes", 500)),
        "gmt_topk_list": gmt_topk_list,
        "gmt_mass_list": gmt_mass_list,
        "gmt_split_signed": bool(_arg(args, "gmt_split_signed", False)),
    }
    if group_name is not None:
        params["group"] = group_name
    return params


def _link(peaks: list[dict[str, object]], genes, args, link_method: str):
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


def _safe_group_name(name: str) -> str:
    out = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
    return out or "group"


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


def _resolve_group_gmt_out_path(group_dir: Path, gmt_out: str | None) -> Path:
    if gmt_out is None or not str(gmt_out).strip():
        return group_dir / "genesets.gmt"
    p = Path(gmt_out)
    if p.is_absolute():
        return group_dir / p.name
    return group_dir / p


def _contrast_peaks(group_vals: list[float], rest_vals: list[float], metric: str, pseudocount: float) -> list[float]:
    out = [0.0] * len(group_vals)
    if metric == "log2fc":
        for i, (x1, x0) in enumerate(zip(group_vals, rest_vals)):
            out[i] = math.log2((float(x1) + pseudocount) / (float(x0) + pseudocount))
        return out
    if metric == "diff":
        for i, (x1, x0) in enumerate(zip(group_vals, rest_vals)):
            out[i] = float(x1) - float(x0)
        return out
    raise ValueError(f"Unsupported contrast_metric: {metric}")


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    peaks, barcodes, matrix_files = read_10x_matrix_dir(args.matrix_dir)
    genes = read_genes_from_gtf(args.gtf)
    groups_tsv = _arg(args, "groups_tsv", None)
    peak_summary = _arg(args, "peak_summary", "sum_counts")
    peak_weight_transform = _arg(args, "peak_weight_transform", "positive")
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
    contrast_metric = _arg(args, "contrast_metric", "log2fc")

    groups = read_groups_tsv(groups_tsv) if groups_tsv else None
    group_indices = make_group_indices(barcodes, groups)
    link_methods = _resolve_link_methods(args)
    primary_link_method = link_methods[0]
    links_by_method = {method: _link(peaks, genes, args, method) for method in link_methods}

    contrast = _resolve_contrast(args)
    pseudocount = _resolve_pseudocount(args)

    files = [
        input_file_record(matrix_files["matrix"], "matrix"),
        input_file_record(matrix_files["barcodes"], "barcodes"),
        input_file_record(matrix_files["peaks_or_features"], "peaks_or_features"),
        input_file_record(args.gtf, "gtf"),
    ]
    if groups_tsv:
        files.append(input_file_record(groups_tsv, "groups_tsv"))

    manifest_rows: list[tuple[str, str, str]] = []
    gene_symbol_by_id = {g.gene_id: g.gene_symbol for g in genes}
    gene_biotype_by_id = {g.gene_id: g.gene_biotype for g in genes}

    if contrast == "group_vs_rest" and not groups_tsv:
        raise ValueError("contrast=group_vs_rest requires --groups_tsv")

    if contrast == "group_vs_rest":
        group_summary, rest_summary, group_sizes = summarize_group_vs_rest(
            Path(matrix_files["matrix"]),
            len(peaks),
            len(barcodes),
            group_indices,
            peak_summary,
        )
    else:
        group_sizes = {g: len(idxs) for g, idxs in group_indices.items()}
        if groups_tsv:
            group_summary = summarize_peaks_by_group(Path(matrix_files["matrix"]), len(peaks), group_indices, peak_summary)
        else:
            group_summary = {
                "all": summarize_peaks(Path(matrix_files["matrix"]), len(peaks), group_indices["all"], peak_summary)
            }
        rest_summary = {}

    unique_output_genes: set[str] = set()
    n_genes_per_group: list[int] = []
    combined_gmt_sets: list[tuple[str, list[str]]] = []
    combined_gmt_plans: list[dict[str, object]] = []

    for group_name, cell_indices in group_indices.items():
        if contrast == "group_vs_rest":
            peak_stat = _contrast_peaks(group_summary[group_name], rest_summary[group_name], contrast_metric, pseudocount)
        else:
            peak_stat = group_summary[group_name]

        full_scores_by_method: dict[str, dict[str, float]] = {}
        for method in link_methods:
            raw_scores = score_genes(peak_stat, links_by_method[method], peak_weight_transform)
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

        full_rows_by_method: dict[str, list[dict[str, object]]] = {}
        for method in link_methods:
            method_rows: list[dict[str, object]] = []
            for rank, gene_id in enumerate(ranked_gene_ids(full_scores_by_method[method]), start=1):
                method_rows.append(
                    {
                        "gene_id": gene_id,
                        "score": float(full_scores_by_method[method][gene_id]),
                        "rank": rank,
                        "gene_symbol": gene_symbol_by_id.get(gene_id),
                        "gene_biotype": gene_biotype_by_id.get(gene_id),
                    }
                )
            full_rows_by_method[method] = method_rows
        full_rows = full_rows_by_method[primary_link_method]

        if groups_tsv:
            group_dir = out_dir / f"group={_safe_group_name(group_name)}"
        else:
            group_dir = out_dir
        group_dir.mkdir(parents=True, exist_ok=True)

        _write_rows(group_dir / "geneset.tsv", selected_rows)
        if emit_full:
            _write_rows(group_dir / "geneset.full.tsv", full_rows)

        group_gmt_path = _resolve_group_gmt_out_path(group_dir, gmt_out)
        group_gmt_sets: list[tuple[str, list[str]]] = []
        group_gmt_plans: list[dict[str, object]] = []
        if emit_gmt:
            for method in link_methods:
                method_sets, method_plans = build_gmt_sets_from_rows(
                    rows=full_rows_by_method[method],
                    base_name=f"atac_sc_10x__group={group_name}__link_method={method}",
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
                    params["link_method"] = method
                    plan["parameters"] = params
                group_gmt_sets.extend(method_sets)
                group_gmt_plans.extend(method_plans)
            write_gmt(group_gmt_sets, group_gmt_path)
            combined_gmt_sets.extend(group_gmt_sets)
            combined_gmt_plans.extend(group_gmt_plans)

        assigned_peaks = len({int(link["peak_index"]) for link in links_by_method[primary_link_method]})
        n_genes_per_group.append(len(selected_rows))
        unique_output_genes.update(gid for gid in selected_gene_ids)

        output_files = [{"path": str(group_dir / "geneset.tsv"), "role": "selected_program"}]
        if emit_full:
            output_files.append({"path": str(group_dir / "geneset.full.tsv"), "role": "full_scores"})
        if emit_gmt:
            output_files.append({"path": str(group_gmt_path), "role": "gmt"})

        contrast_meta: dict[str, object] = {
            "mode": contrast,
            "background_definition": "all_other_cells_in_matrix",
        }
        if contrast == "group_vs_rest":
            contrast_meta.update(
                {
                    "metric": contrast_metric,
                    "pseudocount": pseudocount,
                    "group_size": group_sizes.get(group_name, len(cell_indices)),
                    "rest_size": len(barcodes) - group_sizes.get(group_name, len(cell_indices)),
                }
            )

        meta = make_metadata(
            converter_name="atac_sc_10x",
            parameters=_resolved_parameters(args, group_name),
            data_type="atac_seq",
            assay="single_cell",
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
                "weight_type": "nonnegative" if peak_weight_transform != "signed" else "signed",
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
                "score_definition": "sum over peaks of transformed peak_stat times link_weight",
                "contrast": contrast_meta,
            },
            output_files=output_files,
            gmt={
                "written": emit_gmt,
                "path": (
                    str(group_gmt_path.relative_to(group_dir))
                    if group_gmt_path.is_relative_to(group_dir)
                    else str(group_gmt_path)
                )
                if emit_gmt
                else None,
                "prefer_symbol": gmt_prefer_symbol,
                "require_symbol": gmt_require_symbol,
                "biotype_allowlist": gmt_biotype_allowlist,
                "min_genes": gmt_min_genes,
                "max_genes": gmt_max_genes,
                "plans": group_gmt_plans,
            },
        )
        write_metadata(group_dir / "geneset.meta.json", meta)
        group_rel = str(group_dir.relative_to(out_dir))
        gmt_rel = ""
        if emit_gmt and group_gmt_path.is_relative_to(out_dir):
            gmt_rel = str(group_gmt_path.relative_to(out_dir))
        elif emit_gmt:
            gmt_rel = str(group_gmt_path)
        manifest_rows.append((group_name, group_rel, gmt_rel))

    if groups_tsv:
        root_gmt_path = resolve_gmt_out_path(out_dir, gmt_out)
        if emit_gmt:
            write_gmt(combined_gmt_sets, root_gmt_path)
            root_payload = {
                "groups": len(group_indices),
                "gmt": {
                    "written": True,
                    "path": str(root_gmt_path.relative_to(out_dir))
                    if root_gmt_path.is_relative_to(out_dir)
                    else str(root_gmt_path),
                    "prefer_symbol": gmt_prefer_symbol,
                    "require_symbol": gmt_require_symbol,
                    "biotype_allowlist": gmt_biotype_allowlist,
                    "min_genes": gmt_min_genes,
                    "max_genes": gmt_max_genes,
                    "plans": combined_gmt_plans,
                },
            }
            (out_dir / "manifest.meta.json").write_text(json.dumps(root_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        with (out_dir / "manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["group", "path", "gmt_path"])
            writer.writerows(manifest_rows)

    if groups_tsv:
        return {
            "n_peaks": len(peaks),
            "out_dir": str(out_dir),
            "n_groups": len(group_indices),
            "n_genes_unique": len(unique_output_genes),
            "n_genes_min": min(n_genes_per_group) if n_genes_per_group else 0,
            "n_genes_max": max(n_genes_per_group) if n_genes_per_group else 0,
        }
    return {"n_peaks": len(peaks), "n_genes": n_genes_per_group[0] if n_genes_per_group else 0, "out_dir": str(out_dir), "n_groups": 1}
