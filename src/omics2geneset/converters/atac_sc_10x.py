from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import re
import sys

from omics2geneset.core.gmt import (
    build_gmt_sets_from_rows,
    parse_int_list_csv,
    parse_mass_list_csv,
    parse_str_list_csv,
    resolve_gmt_out_path,
    write_gmt,
)
from omics2geneset.core.atac_programs import (
    CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY,
    CONTRAST_METHOD_NONE,
    PROGRAM_ATLAS_RESIDUAL,
    PROGRAM_DISTAL_ACTIVITY,
    PROGRAM_ENHANCER_BIAS,
    PROGRAM_LINKED_ACTIVITY,
    PROGRAM_PRESET_CONNECTABLE,
    PROGRAM_PROMOTER_ACTIVITY,
    PROGRAM_REF_UBIQUITY_PENALTY,
    PROGRAM_TFIDF_DISTAL,
    compute_peak_idf,
    atlas_residual_scores,
    atlas_stats_for_score_definition,
    contrast_method_enablement_hint,
    default_link_methods_for_preset,
    enhancer_bias_scores,
    link_methods_for_program,
    mask_peak_weights,
    normalize_program_preset,
    peak_values_for_contrast,
    promoter_peak_indices,
    remove_reference_program_methods,
    resolve_auto_contrast_methods,
    resolve_contrast_methods,
    resolve_program_methods,
    score_definition_key,
)
from omics2geneset.core.reference_calibration import peak_overlap_mask, peak_ref_idf_by_overlap
from omics2geneset.core.metadata import input_file_record, make_metadata, write_metadata
from omics2geneset.core.qc import (
    collect_emitted_method_combinations,
    load_marker_genes,
    marker_hit_summary,
    summarize_numeric,
    write_run_summary_files,
)
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
from omics2geneset.io.gtf import read_genes_from_gtf
from omics2geneset.io.mtx_10x import (
    iter_mtx_entries,
    make_group_indices,
    read_cell_metadata_tsv,
    read_10x_matrix_dir,
    read_groups_tsv,
    summarize_condition_within_group,
    summarize_group_vs_rest,
    summarize_peaks,
    summarize_peaks_by_group,
)
from omics2geneset.io.reference_tables import (
    read_atlas_gene_stats_by_score_definition_tsv,
    read_ref_ubiquity_tsv,
)
from omics2geneset.io.region_gene_links import read_region_gene_links_tsv
from omics2geneset.resource_manager import (
    default_resources_dir,
    load_manifest,
    resource_metadata_record,
)


_LINK_METHODS = ("promoter_overlap", "nearest_tss", "distance_decay")
_EXTERNAL_LINK_METHOD = "external"


def _arg(args, name: str, default):
    return getattr(args, name, default)


def _resolve_link_methods(args) -> list[str]:
    link_expr = str(_arg(args, "link_method", "all")).strip()
    preset = normalize_program_preset(_arg(args, "program_preset", PROGRAM_PRESET_CONNECTABLE))
    if not link_expr:
        link_expr = "all"
    if link_expr == "all":
        return [m for m in default_link_methods_for_preset(preset) if m in _LINK_METHODS]
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


def _resolve_contrast(args) -> str:
    contrast = _arg(args, "contrast", None)
    if contrast is not None:
        return str(contrast)
    if normalize_program_preset(_arg(args, "program_preset", PROGRAM_PRESET_CONNECTABLE)) == "connectable":
        if _arg(args, "cell_metadata_tsv", None) and _arg(args, "condition_column", None):
            return "condition_within_group"
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


def _resolve_dataset_label(args) -> str:
    label = _arg(args, "dataset_label", None)
    if label is None:
        return "dataset"
    text = str(label).strip()
    return text or "dataset"


def _resolve_condition_pair(
    args,
    barcodes: list[str],
    metadata_rows: dict[str, dict[str, str]],
) -> tuple[str, str, dict[int, str]]:
    condition_column = str(_arg(args, "condition_column", ""))
    if not condition_column:
        raise ValueError("condition_within_group requires --condition_column")
    condition_by_cell: dict[int, str] = {}
    for idx, barcode in enumerate(barcodes):
        row = metadata_rows.get(barcode)
        if not row:
            continue
        value = row.get(condition_column, "")
        cond = "" if value is None else str(value).strip()
        if not cond:
            continue
        condition_by_cell[idx] = cond

    if not condition_by_cell:
        raise ValueError(
            f"No barcodes in cell_metadata_tsv had non-empty condition values for column: {condition_column}"
        )

    explicit_a = _arg(args, "condition_a", None)
    explicit_b = _arg(args, "condition_b", None)
    unique_conditions = sorted(set(condition_by_cell.values()))

    if explicit_a and explicit_b:
        condition_a = str(explicit_a)
        condition_b = str(explicit_b)
    else:
        if explicit_a:
            condition_a = str(explicit_a)
            others = [x for x in unique_conditions if x != condition_a]
            if len(others) != 1:
                raise ValueError(
                    f"Could not infer --condition_b from condition_column={condition_column}; "
                    "provide both --condition_a and --condition_b"
                )
            condition_b = others[0]
        elif explicit_b:
            condition_b = str(explicit_b)
            others = [x for x in unique_conditions if x != condition_b]
            if len(others) != 1:
                raise ValueError(
                    f"Could not infer --condition_a from condition_column={condition_column}; "
                    "provide both --condition_a and --condition_b"
                )
            condition_a = others[0]
        else:
            if len(unique_conditions) != 2:
                raise ValueError(
                    f"condition_within_group requires exactly two levels in {condition_column} "
                    "(or explicit --condition_a/--condition_b); found: {unique_conditions}"
                )
            condition_a, condition_b = unique_conditions
    if condition_a == condition_b:
        raise ValueError("--condition_a and --condition_b must differ")
    return condition_a, condition_b, condition_by_cell


def _default_ref_ubiquity_resource_id(genome_build: str) -> str | None:
    gb = genome_build.strip().lower()
    if gb in {"hg19", "grch37"}:
        return "ccre_ubiquity_hg19"
    if gb in {"hg38", "grch38"}:
        return "ccre_ubiquity_hg38"
    if gb in {"mm10", "grcm38"}:
        return "ccre_ubiquity_mm10"
    return None


def _default_atlas_resource_id(genome_build: str) -> str | None:
    gb = genome_build.strip().lower()
    if gb in {"hg19", "grch37"}:
        return "atac_reference_profiles_hg19"
    if gb in {"hg38", "grch38"}:
        return "atac_reference_profiles_hg38"
    if gb in {"mm10", "grcm38"}:
        return "atac_reference_profiles_mm10"
    return None


def _resolved_parameters(args, group_name: str | None = None) -> dict[str, object]:
    contrast = _resolve_contrast(args)
    gmt_topk_list = parse_int_list_csv(str(_arg(args, "gmt_topk_list", "200")))
    gmt_mass_list = parse_mass_list_csv(str(_arg(args, "gmt_mass_list", "")))
    gmt_biotype_allowlist = parse_str_list_csv(str(_arg(args, "gmt_biotype_allowlist", "protein_coding")))
    link_methods = _resolve_link_methods(args)
    use_reference_bundle = bool(_arg(args, "use_reference_bundle", True))
    program_preset = normalize_program_preset(_arg(args, "program_preset", PROGRAM_PRESET_CONNECTABLE))
    program_methods = resolve_program_methods(
        "single_cell",
        program_preset,
        _arg(args, "program_methods", None),
    )
    program_methods = remove_reference_program_methods(program_methods)
    contrast_methods = resolve_contrast_methods(
        _arg(args, "contrast_methods", None),
        _arg(args, "program_methods", None),
        use_reference_bundle=use_reference_bundle,
    )
    ref_ubiquity_resource_id = _arg(args, "ref_ubiquity_resource_id", None) or _default_ref_ubiquity_resource_id(
        str(_arg(args, "genome_build", ""))
    )
    atlas_resource_id = _arg(args, "atlas_resource_id", None) or _default_atlas_resource_id(
        str(_arg(args, "genome_build", ""))
    )
    max_distance = (
        _resolve_max_distance(args, link_methods[0])
        if len(link_methods) == 1
        else {method: _resolve_max_distance(args, method) for method in link_methods}
    )
    params: dict[str, object] = {
        "dataset_label": _resolve_dataset_label(args),
        "link_method": _arg(args, "link_method", "all"),
        "link_methods_evaluated": link_methods,
        "primary_link_method": link_methods[0],
        "contrast_methods": _arg(args, "contrast_methods", None),
        "contrast_methods_evaluated": contrast_methods,
        "primary_contrast_method": contrast_methods[0],
        "program_preset": program_preset,
        "program_methods": _arg(args, "program_methods", None),
        "program_methods_evaluated": program_methods,
        "use_reference_bundle": use_reference_bundle,
        "resource_policy": _arg(args, "resource_policy", "skip"),
        "ref_ubiquity_resource_id": ref_ubiquity_resource_id,
        "atlas_resource_id": atlas_resource_id,
        "atlas_metric": _arg(args, "atlas_metric", "zscore"),
        "atlas_eps": _arg(args, "atlas_eps", 1e-6),
        "atlas_min_raw_quantile": _arg(args, "atlas_min_raw_quantile", 0.95),
        "atlas_use_log1p": bool(_arg(args, "atlas_use_log1p", True)),
        "promoter_upstream_bp": _arg(args, "promoter_upstream_bp", 2000),
        "promoter_downstream_bp": _arg(args, "promoter_downstream_bp", 500),
        "max_distance_bp": max_distance,
        "decay_length_bp": _arg(args, "decay_length_bp", 50000),
        "max_genes_per_peak": _arg(args, "max_genes_per_peak", 5),
        "external_linking_enabled": _EXTERNAL_LINK_METHOD in link_methods,
        "peak_summary": _arg(args, "peak_summary", "sum_counts"),
        "peak_weight_transform": _arg(args, "peak_weight_transform", "positive"),
        "contrast": contrast,
        "contrast_metric": _arg(args, "contrast_metric", "log2fc"),
        "contrast_pseudocount": _resolve_pseudocount(args),
        "cell_metadata_tsv_provided": bool(_arg(args, "cell_metadata_tsv", None)),
        "condition_column": _arg(args, "condition_column", None),
        "condition_a": _arg(args, "condition_a", None),
        "condition_b": _arg(args, "condition_b", None),
        "min_cells_per_group": int(_arg(args, "min_cells_per_group", 100)),
        "min_cells_per_condition": int(_arg(args, "min_cells_per_condition", 50)),
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
        "emit_small_gene_sets": bool(_arg(args, "emit_small_gene_sets", False)),
        "marker_qc_enabled": bool(_arg(args, "qc_marker_genes_tsv", None)),
    }
    if group_name is not None:
        params["group"] = group_name
    return params


def _warn_skipped_contrast_methods(
    contrast_methods_skipped: dict[str, str],
    genome_build: str,
    resource_policy: str,
) -> None:
    if not contrast_methods_skipped:
        return
    print("warning: some contrast methods were skipped; continuing with available methods.", file=sys.stderr)
    for method in sorted(contrast_methods_skipped):
        print(
            f"warning: contrast_method={method} skipped: {contrast_methods_skipped[method]}",
            file=sys.stderr,
        )
        hint = contrast_method_enablement_hint(method, genome_build)
        if hint:
            print(f"warning:   {hint}", file=sys.stderr)
    if resource_policy == "skip":
        print(
            "warning: run with --resource_policy fail to require all requested contrast methods.",
            file=sys.stderr,
        )


def _warn_skipped_program_methods(program_methods_skipped: dict[str, str]) -> None:
    if not program_methods_skipped:
        return
    for method in sorted(program_methods_skipped):
        print(
            f"warning: program_method={method} skipped: {program_methods_skipped[method]}",
            file=sys.stderr,
        )


def _warn_gmt_output_diagnostics(gmt_diagnostics: list[dict[str, object]]) -> None:
    for diag in gmt_diagnostics:
        code = str(diag.get("code", ""))
        base_name = str(diag.get("base_name", "unknown"))
        n_genes = diag.get("n_genes", "na")
        reason = str(diag.get("reason", "")).strip()
        suggestion = str(diag.get("suggestion", "")).strip()
        if code in {"small_gene_set_skipped", "no_positive_genes"}:
            print(
                f"warning: skipped GMT output for {base_name}; n_genes={n_genes}. {reason}",
                file=sys.stderr,
            )
        elif code == "small_gene_set_emitted":
            print(
                f"warning: emitted small GMT output for {base_name}; n_genes={n_genes}. {reason}",
                file=sys.stderr,
            )
        elif code == "marginal_gene_count":
            print(
                f"warning: marginal GMT signal for {base_name}; n_genes={n_genes}. {reason}",
                file=sys.stderr,
            )
        else:
            continue
        if suggestion:
            print(f"warning:   {suggestion}", file=sys.stderr)


def _collect_skipped_gmt_outputs(gmt_diagnostics: list[dict[str, object]]) -> list[dict[str, object]]:
    skipped: list[dict[str, object]] = []
    for diag in gmt_diagnostics:
        code = str(diag.get("code", ""))
        if code not in {"small_gene_set_skipped", "no_positive_genes"}:
            continue
        rec: dict[str, object] = {
            "program_method": str(diag.get("program_method", "")),
            "contrast_method": str(diag.get("contrast_method", "")),
            "link_method": str(diag.get("link_method", "")),
            "reason": str(diag.get("reason", "")),
            "n_genes": int(diag.get("n_genes", 0) or 0),
            "code": code,
        }
        direction = str(diag.get("direction", ""))
        if direction:
            rec["direction"] = direction
        skipped.append(rec)
    return skipped


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
    dataset_label = _resolve_dataset_label(args)
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
    gmt_topk_list = parse_int_list_csv(str(_arg(args, "gmt_topk_list", "200")))
    gmt_mass_list = parse_mass_list_csv(str(_arg(args, "gmt_mass_list", "")))
    gmt_split_signed = bool(_arg(args, "gmt_split_signed", False))
    emit_small_gene_sets = bool(_arg(args, "emit_small_gene_sets", False))
    gmt_out = _arg(args, "gmt_out", None)
    marker_genes = load_marker_genes(_arg(args, "qc_marker_genes_tsv", None))
    contrast_metric = _arg(args, "contrast_metric", "log2fc")
    resources_manifest = _arg(args, "resources_manifest", None)
    resources_dir = (
        Path(_arg(args, "resources_dir", "")).expanduser()
        if _arg(args, "resources_dir", None)
        else default_resources_dir()
    )
    resource_policy = str(_arg(args, "resource_policy", "skip"))
    if resource_policy not in {"skip", "fail"}:
        raise ValueError(f"Unsupported resource_policy: {resource_policy}")
    use_reference_bundle = bool(_arg(args, "use_reference_bundle", True))
    program_preset = normalize_program_preset(_arg(args, "program_preset", PROGRAM_PRESET_CONNECTABLE))
    program_methods = resolve_program_methods(
        "single_cell",
        program_preset,
        _arg(args, "program_methods", None),
    )
    program_methods = remove_reference_program_methods(program_methods)
    contrast_methods_requested = resolve_contrast_methods(
        _arg(args, "contrast_methods", None),
        _arg(args, "program_methods", None),
        use_reference_bundle=True,
    )
    contrast_methods = list(contrast_methods_requested)
    program_methods_skipped: dict[str, str] = {}
    contrast_methods_skipped: dict[str, str] = {}
    resources_used: list[dict[str, object]] = []
    if not use_reference_bundle:
        for method in contrast_methods_requested:
            if method == CONTRAST_METHOD_NONE:
                continue
            contrast_methods_skipped.setdefault(
                method,
                "disabled because --use_reference_bundle=false",
            )
        contrast_methods = [m for m in contrast_methods if m == CONTRAST_METHOD_NONE]

    groups = read_groups_tsv(groups_tsv) if groups_tsv else None
    group_indices = make_group_indices(barcodes, groups)
    if not group_indices:
        raise ValueError("No groups were constructed from provided barcodes/groups_tsv")
    link_methods = _resolve_link_methods(args)
    primary_link_method = link_methods[0]

    external_links_unresolved = 0
    if _EXTERNAL_LINK_METHOD in link_methods:
        resolved_external_links, external_links_unresolved = _prepare_external_links(args, genes)
        setattr(args, "_region_gene_links_resolved", resolved_external_links)
    links_by_method = {method: _link(peaks, genes, args, method) for method in link_methods}
    need_promoter = (
        PROGRAM_PROMOTER_ACTIVITY in program_methods
        or PROGRAM_DISTAL_ACTIVITY in program_methods
        or PROGRAM_ENHANCER_BIAS in program_methods
        or PROGRAM_TFIDF_DISTAL in program_methods
    )
    need_distal = (
        PROGRAM_DISTAL_ACTIVITY in program_methods
        or PROGRAM_ENHANCER_BIAS in program_methods
        or PROGRAM_TFIDF_DISTAL in program_methods
    )
    promoter_links = (
        links_by_method["promoter_overlap"]
        if need_promoter and "promoter_overlap" in links_by_method
        else (_link(peaks, genes, args, "promoter_overlap") if need_promoter else [])
    )
    promoter_indices = promoter_peak_indices(promoter_links) if promoter_links else set()
    peak_idf = (
        compute_peak_idf(iter_mtx_entries(Path(matrix_files["matrix"])), len(peaks), len(barcodes))
        if PROGRAM_TFIDF_DISTAL in program_methods
        else []
    )
    needs_ref_ubiquity = (
        PROGRAM_REF_UBIQUITY_PENALTY in contrast_methods
        or CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY in contrast_methods
    )
    needs_atlas = PROGRAM_ATLAS_RESIDUAL in contrast_methods
    manifest_resources: dict[str, dict[str, object]] = {}
    manifest_label = "bundled"
    manifest_warnings: list[str] = []
    ref_peak_idf: list[float] = []
    ref_overlap_count = 0
    atlas_stats_by_definition: dict[str, dict[str, tuple[float, float]]] = {}
    if needs_ref_ubiquity or needs_atlas:
        manifest_label, manifest_resources, _presets, manifest_warnings = load_manifest(resources_manifest)

    if needs_ref_ubiquity:
        ref_resource_id = _arg(args, "ref_ubiquity_resource_id", None) or _default_ref_ubiquity_resource_id(str(args.genome_build))
        try:
            if not ref_resource_id:
                raise ValueError(f"No default ref ubiquity resource id for genome_build={args.genome_build}")
            if ref_resource_id not in manifest_resources:
                raise ValueError(f"Unknown ref ubiquity resource id: {ref_resource_id}")
            ref_entry = manifest_resources[ref_resource_id]
            ref_path = resources_dir / str(ref_entry["filename"])
            if not ref_path.exists():
                raise FileNotFoundError(f"Missing ref ubiquity resource file: {ref_path}")
            ref_rows = read_ref_ubiquity_tsv(ref_path)
            ref_peak_idf = peak_ref_idf_by_overlap(peaks, ref_rows, default_idf=1.0)
            ref_overlap_count = sum(1 for x in peak_overlap_mask(peaks, ref_rows) if x)
            resources_used.append(
                resource_metadata_record(
                    resource_id=ref_resource_id,
                    entry=ref_entry,
                    local_path=ref_path,
                    method=PROGRAM_REF_UBIQUITY_PENALTY,
                )
            )
        except Exception as exc:
            if resource_policy == "fail":
                raise
            contrast_methods_skipped[PROGRAM_REF_UBIQUITY_PENALTY] = str(exc)
            contrast_methods = [m for m in contrast_methods if m != PROGRAM_REF_UBIQUITY_PENALTY]

    if needs_atlas:
        atlas_resource_id = _arg(args, "atlas_resource_id", None) or _default_atlas_resource_id(str(args.genome_build))
        try:
            if not atlas_resource_id:
                raise ValueError(f"No default atlas resource id for genome_build={args.genome_build}")
            if atlas_resource_id not in manifest_resources:
                raise ValueError(f"Unknown atlas resource id: {atlas_resource_id}")
            atlas_entry = manifest_resources[atlas_resource_id]
            atlas_path = resources_dir / str(atlas_entry["filename"])
            if not atlas_path.exists():
                raise FileNotFoundError(f"Missing atlas resource file: {atlas_path}")
            atlas_stats_by_definition = read_atlas_gene_stats_by_score_definition_tsv(atlas_path)
            resources_used.append(
                resource_metadata_record(
                    resource_id=atlas_resource_id,
                    entry=atlas_entry,
                    local_path=atlas_path,
                    method=PROGRAM_ATLAS_RESIDUAL,
                )
            )
        except Exception as exc:
            if resource_policy == "fail":
                raise
            contrast_methods_skipped[PROGRAM_ATLAS_RESIDUAL] = str(exc)
            contrast_methods = [m for m in contrast_methods if m != PROGRAM_ATLAS_RESIDUAL]
    contrast_methods, auto_reason = resolve_auto_contrast_methods(
        contrast_methods,
        ref_ubiquity_ready=bool(ref_peak_idf),
    )
    if auto_reason:
        contrast_methods_skipped[CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY] = auto_reason
    if not contrast_methods:
        contrast_methods = [CONTRAST_METHOD_NONE]
    primary_contrast_method = contrast_methods[0]
    _warn_skipped_contrast_methods(
        contrast_methods_skipped,
        str(args.genome_build),
        resource_policy,
    )

    contrast = _resolve_contrast(args)
    pseudocount = _resolve_pseudocount(args)
    condition_column = str(_arg(args, "condition_column", "")).strip()
    condition_a = ""
    condition_b = ""
    condition_by_cell: dict[int, str] = {}
    condition_counts_by_group: dict[str, dict[str, int]] = {}
    min_cells_per_group = max(1, int(_arg(args, "min_cells_per_group", 100)))
    min_cells_per_condition = max(1, int(_arg(args, "min_cells_per_condition", 50)))

    if contrast == "group_vs_rest" and not groups_tsv:
        raise ValueError("contrast=group_vs_rest requires --groups_tsv")
    if contrast == "condition_within_group":
        cell_metadata_tsv = _arg(args, "cell_metadata_tsv", None)
        explicit_contrast = _arg(args, "contrast", None) is not None
        if not cell_metadata_tsv:
            if explicit_contrast:
                raise ValueError("contrast=condition_within_group requires --cell_metadata_tsv")
            print(
                "warning: condition_within_group contrast requested by preset but no --cell_metadata_tsv "
                "was provided; falling back to group_vs_rest/none.",
                file=sys.stderr,
            )
            contrast = "group_vs_rest" if groups_tsv else "none"
        else:
            metadata_rows = read_cell_metadata_tsv(cell_metadata_tsv)
            try:
                condition_a, condition_b, condition_by_cell = _resolve_condition_pair(args, barcodes, metadata_rows)
            except ValueError as exc:
                if explicit_contrast:
                    raise
                print(
                    "warning: condition_within_group contrast requested by preset but condition labels "
                    f"could not be resolved ({exc}); falling back to group_vs_rest/none.",
                    file=sys.stderr,
                )
                contrast = "group_vs_rest" if groups_tsv else "none"

    files = [
        input_file_record(matrix_files["matrix"], "matrix"),
        input_file_record(matrix_files["barcodes"], "barcodes"),
        input_file_record(matrix_files["peaks_or_features"], "peaks_or_features"),
        input_file_record(args.gtf, "gtf"),
    ]
    if groups_tsv:
        files.append(input_file_record(groups_tsv, "groups_tsv"))
    if contrast == "condition_within_group":
        files.append(input_file_record(_arg(args, "cell_metadata_tsv", None), "cell_metadata_tsv"))
    if _arg(args, "qc_marker_genes_tsv", None):
        files.append(input_file_record(_arg(args, "qc_marker_genes_tsv", None), "qc_marker_genes_tsv"))
    if _EXTERNAL_LINK_METHOD in link_methods:
        files.append(input_file_record(args.region_gene_links_tsv, "region_gene_links_tsv"))
    for r in resources_used:
        files.append(input_file_record(str(r["path"]), f"resource:{r['id']}"))

    manifest_rows: list[tuple[str, str, str]] = []
    gene_symbol_by_id = {g.gene_id: g.gene_symbol for g in genes}
    gene_biotype_by_id = {g.gene_id: g.gene_biotype for g in genes}

    if contrast == "group_vs_rest":
        group_summary, rest_summary, group_sizes = summarize_group_vs_rest(
            Path(matrix_files["matrix"]),
            len(peaks),
            len(barcodes),
            group_indices,
            peak_summary,
        )
        contrast_peak_stats: dict[str, list[float]] = {
            group_name: _contrast_peaks(
                group_summary[group_name],
                rest_summary[group_name],
                contrast_metric,
                pseudocount,
            )
            for group_name in group_indices
        }
        output_groups = []
        for group_name in group_indices:
            group_n = group_sizes.get(group_name, len(group_indices[group_name]))
            if group_n < min_cells_per_group:
                print(
                    "warning: skipping group due to insufficient cells for stable output: "
                    f"group={group_name} n_cells={group_n} "
                    f"(min_cells_per_group={min_cells_per_group})",
                    file=sys.stderr,
                )
                continue
            output_groups.append(group_name)
    else:
        group_sizes = {g: len(idxs) for g, idxs in group_indices.items()}
        if contrast == "condition_within_group":
            condition_summary, condition_counts_by_group = summarize_condition_within_group(
                Path(matrix_files["matrix"]),
                len(peaks),
                group_indices,
                condition_by_cell,
                condition_a,
                condition_b,
                peak_summary,
            )
            contrast_peak_stats = {}
            output_groups = []
            for group_name in group_indices:
                group_n = group_sizes.get(group_name, len(group_indices[group_name]))
                if group_n < min_cells_per_group:
                    print(
                        "warning: skipping group due to insufficient cells for stable output: "
                        f"group={group_name} n_cells={group_n} "
                        f"(min_cells_per_group={min_cells_per_group})",
                        file=sys.stderr,
                    )
                    continue
                counts = condition_counts_by_group[group_name]
                if counts.get(condition_a, 0) < min_cells_per_condition or counts.get(condition_b, 0) < min_cells_per_condition:
                    print(
                        "warning: skipping group because condition labels are too sparse within group: "
                        f"group={group_name} {condition_a}={counts.get(condition_a, 0)} "
                        f"{condition_b}={counts.get(condition_b, 0)} "
                        f"(min_cells_per_condition={min_cells_per_condition})",
                        file=sys.stderr,
                    )
                    continue
                contrast_peak_stats[group_name] = _contrast_peaks(
                    condition_summary[group_name][condition_a],
                    condition_summary[group_name][condition_b],
                    contrast_metric,
                    pseudocount,
                )
                output_groups.append(group_name)
            if not output_groups:
                raise ValueError(
                    "No groups passed min_cells_per_group/min_cells_per_condition for condition_within_group contrast"
                )
        else:
            if groups_tsv:
                group_summary = summarize_peaks_by_group(
                    Path(matrix_files["matrix"]),
                    len(peaks),
                    group_indices,
                    peak_summary,
                )
            else:
                group_summary = {
                    "all": summarize_peaks(
                        Path(matrix_files["matrix"]),
                        len(peaks),
                        group_indices["all"],
                        peak_summary,
                    )
                }
            contrast_peak_stats = {group_name: group_summary[group_name] for group_name in group_indices}
            output_groups = []
            for group_name in group_indices:
                group_n = group_sizes.get(group_name, len(group_indices[group_name]))
                if group_n < min_cells_per_group:
                    print(
                        "warning: skipping group due to insufficient cells for stable output: "
                        f"group={group_name} n_cells={group_n} "
                        f"(min_cells_per_group={min_cells_per_group})",
                        file=sys.stderr,
                    )
                    continue
                output_groups.append(group_name)

    if not output_groups:
        raise ValueError(
            f"No groups passed min_cells_per_group={min_cells_per_group}; "
            "lower --min_cells_per_group to allow small groups."
        )

    unique_output_genes: set[str] = set()
    n_genes_per_group: list[int] = []
    combined_gmt_sets: list[tuple[str, list[str]]] = []
    combined_gmt_plans: list[dict[str, object]] = []
    combined_gmt_diagnostics: list[dict[str, object]] = []

    selected_direction = "OPEN" if peak_weight_transform != "negative" else "CLOSE"
    if contrast == "condition_within_group":
        selected_transform = "positive" if selected_direction == "OPEN" else "negative"
        directions_to_emit = ("OPEN", "CLOSE")
        transform_by_direction = {"OPEN": "positive", "CLOSE": "negative"}
    else:
        selected_direction = "PRIMARY"
        selected_transform = peak_weight_transform
        directions_to_emit = ("PRIMARY",)
        transform_by_direction = {"PRIMARY": peak_weight_transform}

    atlas_metric = str(_arg(args, "atlas_metric", "zscore"))
    atlas_eps = float(_arg(args, "atlas_eps", 1e-6))
    atlas_min_raw_quantile = float(_arg(args, "atlas_min_raw_quantile", 0.95))
    atlas_use_log1p = bool(_arg(args, "atlas_use_log1p", True))
    atlas_missing_score_definitions: set[str] = set()

    for group_name in output_groups:
        cell_indices = group_indices[group_name]
        peak_stat = contrast_peak_stats[group_name]
        direction_scores_by_contrast_by_method: dict[str, dict[str, dict[str, dict[str, float]]]] = {}
        full_rows_by_contrast_by_method_by_direction: dict[str, dict[str, dict[str, list[dict[str, object]]]]] = {}
        additional_program_rows_by_contrast_by_direction: dict[
            str, dict[str, dict[str, dict[str, list[dict[str, object]]]]]
        ] = {}

        for direction in directions_to_emit:
            transform_mode = transform_by_direction[direction]
            def _linked_scores_for_peak_values(cur_peak_values: list[float]) -> dict[str, dict[str, float]]:
                out: dict[str, dict[str, float]] = {}
                for method in link_methods:
                    raw_scores = score_genes(cur_peak_values, links_by_method[method], transform_mode)
                    out[method] = {g: float(s) for g, s in raw_scores.items() if float(s) != 0.0}
                return out

            def _family_scores_for_peak_values(cur_peak_values: list[float]) -> dict[str, dict[str, dict[str, float]]]:
                out: dict[str, dict[str, dict[str, float]]] = {}
                for method in link_methods:
                    out[method] = {}
                    promoter_scores: dict[str, float] = {}
                    distal_scores: dict[str, float] = {}
                    if PROGRAM_PROMOTER_ACTIVITY in program_methods or PROGRAM_ENHANCER_BIAS in program_methods:
                        promoter_peak_stat = mask_peak_weights(cur_peak_values, include_indices=promoter_indices)
                        promoter_raw = score_genes(promoter_peak_stat, links_by_method[method], transform_mode)
                        promoter_scores = {g: float(s) for g, s in promoter_raw.items() if float(s) != 0.0}
                        if PROGRAM_PROMOTER_ACTIVITY in program_methods:
                            out[method][PROGRAM_PROMOTER_ACTIVITY] = promoter_scores
                    if PROGRAM_DISTAL_ACTIVITY in program_methods or PROGRAM_ENHANCER_BIAS in program_methods:
                        distal_peak_stat = mask_peak_weights(cur_peak_values, exclude_indices=promoter_indices)
                        distal_raw = score_genes(distal_peak_stat, links_by_method[method], transform_mode)
                        distal_scores = {g: float(s) for g, s in distal_raw.items() if float(s) != 0.0}
                        if PROGRAM_DISTAL_ACTIVITY in program_methods:
                            out[method][PROGRAM_DISTAL_ACTIVITY] = distal_scores
                    if PROGRAM_ENHANCER_BIAS in program_methods:
                        bias_scores = enhancer_bias_scores(promoter_scores, distal_scores)
                        out[method][PROGRAM_ENHANCER_BIAS] = {g: float(s) for g, s in bias_scores.items() if float(s) != 0.0}
                    if PROGRAM_TFIDF_DISTAL in program_methods:
                        tfidf_peak_stat = [float(cur_peak_values[i]) * float(peak_idf[i]) for i in range(len(cur_peak_values))]
                        tfidf_distal_peak_stat = mask_peak_weights(tfidf_peak_stat, exclude_indices=promoter_indices)
                        tfidf_distal_raw = score_genes(tfidf_distal_peak_stat, links_by_method[method], transform_mode)
                        out[method][PROGRAM_TFIDF_DISTAL] = {
                            g: float(s) for g, s in tfidf_distal_raw.items() if float(s) != 0.0
                        }
                return out

            full_scores_none_by_method = _linked_scores_for_peak_values(peak_stat)
            family_scores_none_by_method = _family_scores_for_peak_values(peak_stat)
            full_scores_by_contrast_by_method: dict[str, dict[str, dict[str, float]]] = {}
            family_scores_by_contrast_by_method: dict[str, dict[str, dict[str, dict[str, float]]]] = {}
            for contrast_method in contrast_methods:
                if contrast_method == PROGRAM_ATLAS_RESIDUAL:
                    full_scores_by_contrast_by_method[contrast_method] = {}
                    family_scores_by_contrast_by_method[contrast_method] = {}
                    for method in link_methods:
                        linked_key = score_definition_key(PROGRAM_LINKED_ACTIVITY, method)
                        linked_stats = atlas_stats_for_score_definition(atlas_stats_by_definition, linked_key)
                        if linked_stats is None:
                            if resource_policy == "fail":
                                raise ValueError(
                                    "Missing atlas baseline for score_definition="
                                    f"{linked_key}. Provide matching score definitions in atlas resource."
                                )
                            atlas_missing_score_definitions.add(linked_key)
                            full_scores_by_contrast_by_method[contrast_method][method] = {}
                        else:
                            atlas_scores = atlas_residual_scores(
                                full_scores_none_by_method[method],
                                linked_stats,
                                atlas_metric,
                                atlas_eps,
                                min_raw_quantile=atlas_min_raw_quantile,
                                use_log1p=atlas_use_log1p,
                            )
                            full_scores_by_contrast_by_method[contrast_method][method] = {
                                g: float(s) for g, s in atlas_scores.items() if float(s) != 0.0
                            }
                        family_scores_by_contrast_by_method[contrast_method][method] = {}
                        for program_method, base_scores in family_scores_none_by_method.get(method, {}).items():
                            score_key = score_definition_key(program_method, method)
                            score_stats = atlas_stats_for_score_definition(atlas_stats_by_definition, score_key)
                            if score_stats is None:
                                if resource_policy == "fail":
                                    raise ValueError(
                                        "Missing atlas baseline for score_definition="
                                        f"{score_key}. Provide matching score definitions in atlas resource."
                                    )
                                atlas_missing_score_definitions.add(score_key)
                                family_scores_by_contrast_by_method[contrast_method][method][program_method] = {}
                                continue
                            residual_scores = atlas_residual_scores(
                                base_scores,
                                score_stats,
                                atlas_metric,
                                atlas_eps,
                                min_raw_quantile=atlas_min_raw_quantile,
                                use_log1p=atlas_use_log1p,
                            )
                            family_scores_by_contrast_by_method[contrast_method][method][program_method] = {
                                g: float(s) for g, s in residual_scores.items() if float(s) != 0.0
                            }
                    continue

                contrast_peak_values = peak_values_for_contrast(
                    peak_stat,
                    contrast_method,
                    ref_peak_idf if ref_peak_idf else None,
                )
                full_scores_by_contrast_by_method[contrast_method] = _linked_scores_for_peak_values(contrast_peak_values)
                family_scores_by_contrast_by_method[contrast_method] = _family_scores_for_peak_values(contrast_peak_values)

            direction_scores_by_contrast_by_method[direction] = full_scores_by_contrast_by_method

            rows_for_direction: dict[str, dict[str, list[dict[str, object]]]] = {}
            for contrast_method in contrast_methods:
                rows_for_direction[contrast_method] = {}
                for method in link_methods:
                    rows_for_direction[contrast_method][method] = _rows_from_scores(
                        full_scores_by_contrast_by_method[contrast_method][method],
                        gene_symbol_by_id,
                        gene_biotype_by_id,
                    )
            full_rows_by_contrast_by_method_by_direction[direction] = rows_for_direction

            program_rows_by_contrast_by_method: dict[
                str, dict[str, dict[str, list[dict[str, object]]]]
            ] = {}
            for contrast_method in contrast_methods:
                rows_by_method: dict[str, dict[str, list[dict[str, object]]]] = {}
                for method in link_methods:
                    rows_by_program: dict[str, list[dict[str, object]]] = {}
                    for program_method in program_methods:
                        if program_method == PROGRAM_LINKED_ACTIVITY:
                            continue
                        rows_by_program[program_method] = _rows_from_scores(
                            family_scores_by_contrast_by_method.get(contrast_method, {})
                            .get(method, {})
                            .get(program_method, {}),
                            gene_symbol_by_id,
                            gene_biotype_by_id,
                        )
                    rows_by_method[method] = rows_by_program
                program_rows_by_contrast_by_method[contrast_method] = rows_by_method
            additional_program_rows_by_contrast_by_direction[direction] = program_rows_by_contrast_by_method

        full_scores = direction_scores_by_contrast_by_method[
            selected_direction if selected_direction in direction_scores_by_contrast_by_method else "PRIMARY"
        ][primary_contrast_method][primary_link_method]
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

        full_rows = full_rows_by_contrast_by_method_by_direction[
            selected_direction if selected_direction in full_rows_by_contrast_by_method_by_direction else "PRIMARY"
        ][primary_contrast_method][primary_link_method]

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
        group_requested_gmt_outputs: list[dict[str, str]] = []
        group_gmt_diagnostics: list[dict[str, object]] = []
        linked_output_methods = link_methods_for_program(program_preset, PROGRAM_LINKED_ACTIVITY, link_methods)
        if PROGRAM_LINKED_ACTIVITY in program_methods and not linked_output_methods:
            program_methods_skipped[PROGRAM_LINKED_ACTIVITY] = (
                "no compatible link_method selected for preset; expected nearest_tss for connectable/default"
            )
        if emit_gmt:
            for direction in directions_to_emit:
                direction_suffix = (
                    f"__direction={direction}"
                    if contrast == "condition_within_group"
                    else ""
                )
                contrast_suffix = f"__contrast={contrast}"
                if contrast == "condition_within_group":
                    contrast_suffix += (
                        f"__condition_column={condition_column}"
                        f"__condition_a={condition_a}"
                        f"__condition_b={condition_b}"
                    )
                base_prefix = (
                    f"atac_sc_10x__dataset={dataset_label}__group={group_name}"
                    f"{contrast_suffix}{direction_suffix}"
                )

                if PROGRAM_LINKED_ACTIVITY in program_methods:
                    for contrast_method in contrast_methods:
                        for method in linked_output_methods:
                            group_requested_gmt_outputs.append(
                                {
                                    "program_method": PROGRAM_LINKED_ACTIVITY,
                                    "contrast_method": contrast_method,
                                    "link_method": method,
                                    "direction": direction,
                                }
                            )
                            method_sets, method_plans = build_gmt_sets_from_rows(
                                rows=full_rows_by_contrast_by_method_by_direction[direction][contrast_method][method],
                                base_name=(
                                    f"{base_prefix}__program={PROGRAM_LINKED_ACTIVITY}"
                                    f"__contrast_method={contrast_method}"
                                    f"__link_method={method}"
                                ),
                                prefer_symbol=gmt_prefer_symbol,
                                min_genes=gmt_min_genes,
                                max_genes=gmt_max_genes,
                                topk_list=gmt_topk_list,
                                mass_list=gmt_mass_list,
                                split_signed=gmt_split_signed,
                                require_symbol=gmt_require_symbol,
                                allowed_biotypes={b.lower() for b in gmt_biotype_allowlist} if gmt_biotype_allowlist else None,
                                emit_small_gene_sets=emit_small_gene_sets,
                                diagnostics=group_gmt_diagnostics,
                                context={
                                    "group": group_name,
                                    "program_method": PROGRAM_LINKED_ACTIVITY,
                                    "contrast_method": contrast_method,
                                    "link_method": method,
                                    "direction": direction,
                                },
                            )
                            for plan in method_plans:
                                params = dict(plan.get("parameters", {}))
                                params["program_method"] = PROGRAM_LINKED_ACTIVITY
                                params["contrast_method"] = contrast_method
                                params["link_method"] = method
                                params["direction"] = direction
                                plan["parameters"] = params
                            group_gmt_sets.extend(method_sets)
                            group_gmt_plans.extend(method_plans)

                for contrast_method, rows_by_method in additional_program_rows_by_contrast_by_direction[direction].items():
                    for method, rows_by_program in rows_by_method.items():
                        for program_method, rows in rows_by_program.items():
                            allowed_methods = link_methods_for_program(program_preset, program_method, link_methods)
                            if method not in allowed_methods:
                                continue
                            group_requested_gmt_outputs.append(
                                {
                                    "program_method": program_method,
                                    "contrast_method": contrast_method,
                                    "link_method": method,
                                    "direction": direction,
                                }
                            )
                            method_sets, method_plans = build_gmt_sets_from_rows(
                                rows=rows,
                                base_name=(
                                    f"{base_prefix}__program={program_method}"
                                    f"__contrast_method={contrast_method}"
                                    f"__link_method={method}"
                                ),
                                prefer_symbol=gmt_prefer_symbol,
                                min_genes=gmt_min_genes,
                                max_genes=gmt_max_genes,
                                topk_list=gmt_topk_list,
                                mass_list=gmt_mass_list,
                                split_signed=gmt_split_signed,
                                require_symbol=gmt_require_symbol,
                                allowed_biotypes={b.lower() for b in gmt_biotype_allowlist} if gmt_biotype_allowlist else None,
                                emit_small_gene_sets=emit_small_gene_sets,
                                diagnostics=group_gmt_diagnostics,
                                context={
                                    "group": group_name,
                                    "program_method": program_method,
                                    "contrast_method": contrast_method,
                                    "link_method": method,
                                    "direction": direction,
                                },
                            )
                            for plan in method_plans:
                                params = dict(plan.get("parameters", {}))
                                params["program_method"] = program_method
                                params["contrast_method"] = contrast_method
                                params["link_method"] = method
                                params["direction"] = direction
                                plan["parameters"] = params
                            group_gmt_sets.extend(method_sets)
                            group_gmt_plans.extend(method_plans)
            for program_method in program_methods:
                if program_method == PROGRAM_LINKED_ACTIVITY:
                    continue
                allowed_methods = link_methods_for_program(program_preset, program_method, link_methods)
                if allowed_methods:
                    continue
                if program_method not in program_methods_skipped:
                    program_methods_skipped[program_method] = "no compatible link_method selected for active preset"
            write_gmt(group_gmt_sets, group_gmt_path)
            combined_gmt_sets.extend(group_gmt_sets)
            combined_gmt_plans.extend(group_gmt_plans)
            combined_gmt_diagnostics.extend(group_gmt_diagnostics)

        _warn_gmt_output_diagnostics(group_gmt_diagnostics)

        assigned_peaks = len({int(link["peak_index"]) for link in links_by_method[primary_link_method]})
        n_genes_per_group.append(len(selected_rows))
        unique_output_genes.update(gid for gid in selected_gene_ids)
        link_assignment: dict[str, dict[str, float]] = {}
        for method in link_methods:
            assigned = len({int(link["peak_index"]) for link in links_by_method[method]})
            link_assignment[method] = {
                "n_features_assigned": int(assigned),
                "fraction_features_assigned": float(assigned / len(peaks) if peaks else 0.0),
            }
        marker_qc = marker_hit_summary(selected_rows, marker_genes)
        emitted_combinations = collect_emitted_method_combinations(group_gmt_plans)
        skipped_gmt_outputs = _collect_skipped_gmt_outputs(group_gmt_diagnostics)
        if not emitted_combinations:
            fallback_direction = selected_direction if contrast == "condition_within_group" else "PRIMARY"
            emitted_combinations = [
                {
                    "program_method": PROGRAM_LINKED_ACTIVITY,
                    "contrast_method": primary_contrast_method,
                    "link_method": primary_link_method,
                    "direction": fallback_direction,
                }
            ]
        ref_summary: dict[str, object] | None = None
        if ref_peak_idf:
            adjusted = [float(w) * float(i) for w, i in zip(peak_stat, ref_peak_idf)]
            ref_summary = {
                "n_overlapping_peaks": ref_overlap_count,
                "fraction_overlapping_peaks": (ref_overlap_count / len(peaks) if peaks else 0.0),
                "idf_stats": summarize_numeric([float(x) for x in ref_peak_idf]),
                "adjusted_peak_weight_stats": summarize_numeric(adjusted),
            }

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
        elif contrast == "condition_within_group":
            contrast_meta.update(
                {
                    "contrast_type": "condition_within_group",
                    "condition_column": condition_column,
                    "condition_a": condition_a,
                    "condition_b": condition_b,
                    "metric": contrast_metric,
                    "pseudocount": pseudocount,
                    "n_cells": len(cell_indices),
                    "n_cells_a": condition_counts_by_group.get(group_name, {}).get(condition_a, 0),
                    "n_cells_b": condition_counts_by_group.get(group_name, {}).get(condition_b, 0),
                    "directions_emitted": ["OPEN", "CLOSE"],
                    "selected_direction": selected_direction,
                }
            )

        run_summary_payload: dict[str, object] = {
            "converter": "atac_sc_10x",
            "dataset_label": dataset_label,
            "group": group_name,
            "program_preset": program_preset,
            "primary_program_method": PROGRAM_LINKED_ACTIVITY,
            "primary_link_method": primary_link_method,
            "primary_contrast_method": primary_contrast_method,
            "selected_direction": selected_direction if contrast == "condition_within_group" else "PRIMARY",
            "n_input_peaks": len(peaks),
            "n_cells": len(cell_indices),
            "link_assignment": link_assignment,
            "promoter_peak_count": len(promoter_indices),
            "distal_peak_count": len(peaks) - len(promoter_indices),
            "emitted_method_combinations": emitted_combinations,
            "requested_gmt_outputs": group_requested_gmt_outputs,
            "skipped_gmt_outputs": skipped_gmt_outputs,
            "program_methods_skipped": program_methods_skipped,
            "contrast_methods_skipped": contrast_methods_skipped,
        }
        if marker_qc is not None:
            run_summary_payload["marker_qc"] = marker_qc
        if ref_summary is not None:
            run_summary_payload["ref_ubiquity"] = ref_summary
        run_summary_json_path, run_summary_txt_path = write_run_summary_files(group_dir, run_summary_payload)
        output_files.append({"path": str(run_summary_json_path), "role": "run_summary_json"})
        output_files.append({"path": str(run_summary_txt_path), "role": "run_summary_text"})

        params = _resolved_parameters(args, group_name)
        params["program_methods_active"] = program_methods
        params["program_methods_skipped"] = program_methods_skipped
        params["contrast_methods_active"] = contrast_methods
        params["contrast_methods_skipped"] = contrast_methods_skipped
        params["primary_contrast_method"] = primary_contrast_method
        params["dataset_label"] = dataset_label
        if contrast == "condition_within_group":
            params["condition_column"] = condition_column
            params["condition_a"] = condition_a
            params["condition_b"] = condition_b

        summary_payload: dict[str, object] = {
            "n_input_features": len(peaks),
            "n_genes": len(selected_rows),
            "n_features_assigned": assigned_peaks,
            "fraction_features_assigned": assigned_peaks / len(peaks) if peaks else 0.0,
            "n_external_links_unresolved_gene_id": external_links_unresolved,
            "n_program_methods": len(program_methods),
            "n_contrast_methods": len(contrast_methods),
            "n_contrast_methods_skipped": len(contrast_methods_skipped),
            "n_resource_manifest_warnings": len(manifest_warnings),
        }
        if marker_qc is not None:
            summary_payload["marker_qc"] = marker_qc

        meta = make_metadata(
            converter_name="atac_sc_10x",
            parameters=params,
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
            summary=summary_payload,
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
                "contrast_methods": contrast_methods,
                "contrast_methods_skipped": contrast_methods_skipped,
                "primary_contrast_method": primary_contrast_method,
                "program_methods": program_methods,
                "program_methods_skipped": program_methods_skipped,
                "primary_program_method": PROGRAM_LINKED_ACTIVITY,
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
                "emit_small_gene_sets": emit_small_gene_sets,
                "requested_outputs": group_requested_gmt_outputs,
                "emitted_outputs": emitted_combinations,
                "skipped_outputs": skipped_gmt_outputs,
                "diagnostics": group_gmt_diagnostics,
                "plans": group_gmt_plans,
            },
        )
        if resources_used or manifest_warnings:
            meta["resources"] = {
                "manifest": manifest_label,
                "resources_dir": str(resources_dir),
                "used": resources_used,
                "warnings": manifest_warnings,
            }
        write_metadata(group_dir / "geneset.meta.json", meta)
        group_rel = str(group_dir.relative_to(out_dir))
        gmt_rel = ""
        if emit_gmt and group_gmt_path.is_relative_to(out_dir):
            gmt_rel = str(group_gmt_path.relative_to(out_dir))
        elif emit_gmt:
            gmt_rel = str(group_gmt_path)
        manifest_rows.append((group_name, group_rel, gmt_rel))

    _warn_skipped_program_methods(program_methods_skipped)

    if atlas_missing_score_definitions:
        msg = (
            "atlas_residual baseline missing for score definitions: "
            + ", ".join(sorted(atlas_missing_score_definitions))
        )
        contrast_methods_skipped.setdefault(PROGRAM_ATLAS_RESIDUAL, msg)
        print(f"warning: {msg}", file=sys.stderr)
        print(
            "warning: provide per-score-definition atlas baselines or run with --contrast_methods none/ref_ubiquity_penalty.",
            file=sys.stderr,
        )

    if groups_tsv:
        root_gmt_path = resolve_gmt_out_path(out_dir, gmt_out)
        if emit_gmt:
            write_gmt(combined_gmt_sets, root_gmt_path)
            root_payload = {
                "groups": len(output_groups),
                "program_methods_skipped": program_methods_skipped,
                "contrast_methods_skipped": contrast_methods_skipped,
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
                    "emit_small_gene_sets": emit_small_gene_sets,
                    "diagnostics": combined_gmt_diagnostics,
                    "plans": combined_gmt_plans,
                },
            }
            if resources_used or manifest_warnings:
                root_payload["resources"] = {
                    "manifest": manifest_label,
                    "resources_dir": str(resources_dir),
                    "used": resources_used,
                    "warnings": manifest_warnings,
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
            "n_groups": len(output_groups),
            "n_genes_unique": len(unique_output_genes),
            "n_genes_min": min(n_genes_per_group) if n_genes_per_group else 0,
            "n_genes_max": max(n_genes_per_group) if n_genes_per_group else 0,
            "program_methods": program_methods,
            "program_methods_skipped": program_methods_skipped,
            "contrast_methods": contrast_methods,
            "contrast_methods_skipped": contrast_methods_skipped,
        }
    return {
        "n_peaks": len(peaks),
        "n_genes": n_genes_per_group[0] if n_genes_per_group else 0,
        "out_dir": str(out_dir),
        "n_groups": 1,
        "program_methods": program_methods,
        "program_methods_skipped": program_methods_skipped,
        "contrast_methods": contrast_methods,
        "contrast_methods_skipped": contrast_methods_skipped,
    }
