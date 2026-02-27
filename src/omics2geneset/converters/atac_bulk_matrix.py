from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

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
    atlas_residual_scores,
    atlas_stats_for_score_definition,
    calibration_method_enablement_hint,
    default_link_methods_for_preset,
    enhancer_bias_scores,
    link_methods_for_program,
    mask_peak_weights,
    normalize_program_preset,
    peak_values_for_contrast,
    promoter_peak_indices,
    remove_reference_program_methods,
    resolve_auto_calibration_methods,
    resolve_calibration_methods,
    resolve_program_methods,
    score_definition_key,
)
from omics2geneset.core.gmt import (
    build_gmt_sets_from_rows,
    parse_int_list_csv,
    parse_mass_list_csv,
    parse_str_list_csv,
    resolve_gmt_out_path,
    write_gmt,
)
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
from omics2geneset.core.reference_calibration import peak_overlap_mask, peak_ref_idf_by_overlap
from omics2geneset.core.scoring import score_genes
from omics2geneset.core.selection import (
    global_l1_weights,
    ranked_gene_ids,
    select_quantile,
    select_threshold,
    select_top_k,
    within_set_l1_weights,
)
from omics2geneset.io.bed import read_bed
from omics2geneset.io.gtf import read_genes_from_gtf
from omics2geneset.io.reference_tables import (
    read_atlas_gene_stats_by_score_definition_tsv,
    read_ref_ubiquity_tsv,
)
from omics2geneset.io.region_gene_links import read_region_gene_links_tsv
from omics2geneset.resource_manager import default_resources_dir, load_manifest, resource_metadata_record


_LINK_METHODS = ("promoter_overlap", "nearest_tss", "distance_decay")
_EXTERNAL_LINK_METHOD = "external"


def _arg(args, name: str, default):
    aliases = {
        "calibration_methods": ("contrast_methods",),
        "contrast_methods": ("calibration_methods",),
        "study_contrast": ("contrast",),
        "contrast": ("study_contrast",),
    }
    if hasattr(args, name):
        value = getattr(args, name)
        if value is not None:
            return value
    for alias in aliases.get(name, ()):
        if hasattr(args, alias):
            value = getattr(args, alias)
            if value is not None:
                return value
    return default


def _resolve_dataset_label(args) -> str:
    label = _arg(args, "dataset_label", None)
    if label is None:
        return "dataset"
    text = str(label).strip()
    return text or "dataset"


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


def _resolved_parameters(
    args,
    condition_a: str,
    condition_b: str,
) -> dict[str, object]:
    gmt_topk_list = parse_int_list_csv(str(_arg(args, "gmt_topk_list", "200")))
    gmt_mass_list = parse_mass_list_csv(str(_arg(args, "gmt_mass_list", "")))
    gmt_biotype_allowlist = parse_str_list_csv(str(_arg(args, "gmt_biotype_allowlist", "protein_coding")))
    link_methods = _resolve_link_methods(args)
    use_reference_bundle = bool(_arg(args, "use_reference_bundle", True))
    program_preset = normalize_program_preset(_arg(args, "program_preset", PROGRAM_PRESET_CONNECTABLE))
    program_methods = resolve_program_methods(
        "bulk",
        program_preset,
        _arg(args, "program_methods", None),
    )
    program_methods = remove_reference_program_methods(program_methods)
    calibration_methods = resolve_calibration_methods(
        _arg(args, "calibration_methods", None),
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
    return {
        "dataset_label": _resolve_dataset_label(args),
        "link_method": _arg(args, "link_method", "all"),
        "link_methods_evaluated": link_methods,
        "primary_link_method": link_methods[0],
        "study_contrast": "condition_between_samples",
        "calibration_methods": _arg(args, "calibration_methods", None),
        "calibration_methods_evaluated": calibration_methods,
        "primary_calibration_method": calibration_methods[0],
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
        "study_contrast": "condition_between_samples",
        "contrast": "condition_between_samples",
        "contrast_metric": _arg(args, "contrast_metric", "log2fc"),
        "contrast_pseudocount": float(_arg(args, "contrast_pseudocount", 1.0)),
        "sample_id_column": _arg(args, "sample_id_column", "sample_id"),
        "condition_column": _arg(args, "condition_column", "condition"),
        "condition_a": condition_a,
        "condition_b": condition_b,
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


def _warn_skipped_calibration_methods(
    calibration_methods_skipped: dict[str, str],
    genome_build: str,
    resource_policy: str,
) -> None:
    if not calibration_methods_skipped:
        return
    print("warning: some calibration methods were skipped; continuing with available methods.", file=sys.stderr)
    for method in sorted(calibration_methods_skipped):
        print(
            f"warning: calibration_method={method} skipped: {calibration_methods_skipped[method]}",
            file=sys.stderr,
        )
        hint = calibration_method_enablement_hint(method, genome_build)
        if hint:
            print(f"warning:   {hint}", file=sys.stderr)
    if resource_policy == "skip":
        print(
            "warning: run with --resource_policy fail to require all requested calibration methods.",
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
            "calibration_method": str(diag.get("calibration_method", "")),
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


def _read_peak_matrix_tsv(path: str | Path) -> tuple[list[str], list[list[float]]]:
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("peak_matrix_tsv is missing header")
        reserved = {"chrom", "start", "end", "peak_id", "name"}
        sample_ids = [h for h in reader.fieldnames if h and h not in reserved]
        if not sample_ids:
            raise ValueError("peak_matrix_tsv must include one or more sample columns")
        rows: list[list[float]] = []
        for row in reader:
            rows.append([float(row[sample]) for sample in sample_ids])
    return sample_ids, rows


def _read_sample_conditions(
    path: str | Path,
    sample_id_column: str,
    condition_column: str,
) -> dict[str, str]:
    out: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("sample_metadata_tsv is missing header")
        fields = set(reader.fieldnames)
        if sample_id_column not in fields:
            raise ValueError(f"sample_metadata_tsv missing required column: {sample_id_column}")
        if condition_column not in fields:
            raise ValueError(f"sample_metadata_tsv missing required column: {condition_column}")
        for row in reader:
            sample_id = str(row.get(sample_id_column, "")).strip()
            condition = str(row.get(condition_column, "")).strip()
            if not sample_id or not condition:
                continue
            out[sample_id] = condition
    return out


def _resolve_condition_pair(
    args,
    sample_ids: list[str],
    sample_conditions: dict[str, str],
) -> tuple[str, str, list[int], list[int]]:
    samples_with_condition = [(idx, sid, sample_conditions.get(sid, "")) for idx, sid in enumerate(sample_ids)]
    unique_conditions = sorted({cond for _idx, _sid, cond in samples_with_condition if cond})
    if not unique_conditions:
        raise ValueError("No condition labels matched sample IDs in peak_matrix_tsv")

    explicit_a = _arg(args, "condition_a", None)
    explicit_b = _arg(args, "condition_b", None)
    if explicit_a and explicit_b:
        condition_a = str(explicit_a)
        condition_b = str(explicit_b)
    else:
        if explicit_a:
            condition_a = str(explicit_a)
            others = [x for x in unique_conditions if x != condition_a]
            if len(others) != 1:
                raise ValueError("Could not infer --condition_b; provide both --condition_a and --condition_b")
            condition_b = others[0]
        elif explicit_b:
            condition_b = str(explicit_b)
            others = [x for x in unique_conditions if x != condition_b]
            if len(others) != 1:
                raise ValueError("Could not infer --condition_a; provide both --condition_a and --condition_b")
            condition_a = others[0]
        else:
            if len(unique_conditions) != 2:
                raise ValueError(
                    "Need exactly two condition levels in sample_metadata_tsv, "
                    "or pass --condition_a and --condition_b"
                )
            condition_a, condition_b = unique_conditions
    if condition_a == condition_b:
        raise ValueError("--condition_a and --condition_b must differ")

    a_indices = [idx for idx, _sid, cond in samples_with_condition if cond == condition_a]
    b_indices = [idx for idx, _sid, cond in samples_with_condition if cond == condition_b]
    if not a_indices or not b_indices:
        raise ValueError(
            f"Resolved conditions had empty groups: {condition_a}={len(a_indices)} {condition_b}={len(b_indices)}"
        )
    return condition_a, condition_b, a_indices, b_indices


def _compute_peak_contrast(
    matrix_rows: list[list[float]],
    a_indices: list[int],
    b_indices: list[int],
    metric: str,
    pseudocount: float,
) -> list[float]:
    stats: list[float] = []
    n_a = float(len(a_indices))
    n_b = float(len(b_indices))
    for row in matrix_rows:
        mean_a = sum(float(row[i]) for i in a_indices) / n_a
        mean_b = sum(float(row[i]) for i in b_indices) / n_b
        if metric == "log2fc":
            stats.append(math.log2((mean_a + pseudocount) / (mean_b + pseudocount)))
        elif metric == "diff":
            stats.append(mean_a - mean_b)
        else:
            raise ValueError(f"Unsupported contrast_metric: {metric}")
    return stats


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    peaks = read_bed(args.peak_bed)
    genes = read_genes_from_gtf(args.gtf)
    sample_ids, matrix_rows = _read_peak_matrix_tsv(args.peak_matrix_tsv)
    if len(matrix_rows) != len(peaks):
        raise ValueError(
            "peak_matrix_tsv row count must match peak_bed row count: "
            f"{len(matrix_rows)} vs {len(peaks)}"
        )

    sample_id_column = str(_arg(args, "sample_id_column", "sample_id"))
    condition_column = str(_arg(args, "condition_column", "condition"))
    sample_conditions = _read_sample_conditions(args.sample_metadata_tsv, sample_id_column, condition_column)
    condition_a, condition_b, a_indices, b_indices = _resolve_condition_pair(args, sample_ids, sample_conditions)

    contrast_metric = str(_arg(args, "contrast_metric", "log2fc"))
    pseudocount = float(_arg(args, "contrast_pseudocount", 1.0))
    peak_stat = _compute_peak_contrast(matrix_rows, a_indices, b_indices, contrast_metric, pseudocount)

    peak_weight_transform = str(_arg(args, "peak_weight_transform", "positive"))
    selected_direction = "OPEN" if peak_weight_transform != "negative" else "CLOSE"
    selected_transform = "positive" if selected_direction == "OPEN" else "negative"
    marker_genes = load_marker_genes(_arg(args, "qc_marker_genes_tsv", None))

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
    dataset_label = _resolve_dataset_label(args)

    resource_policy = str(_arg(args, "resource_policy", "skip"))
    if resource_policy not in {"skip", "fail"}:
        raise ValueError(f"Unsupported resource_policy: {resource_policy}")
    use_reference_bundle = bool(_arg(args, "use_reference_bundle", True))
    program_preset = normalize_program_preset(_arg(args, "program_preset", PROGRAM_PRESET_CONNECTABLE))
    program_methods = resolve_program_methods(
        "bulk",
        program_preset,
        _arg(args, "program_methods", None),
    )
    program_methods = remove_reference_program_methods(program_methods)
    calibration_methods_requested = resolve_calibration_methods(
        _arg(args, "calibration_methods", None),
        _arg(args, "program_methods", None),
        use_reference_bundle=True,
    )
    calibration_methods = list(calibration_methods_requested)
    program_methods_skipped: dict[str, str] = {}
    calibration_methods_skipped: dict[str, str] = {}
    resources_used: list[dict[str, object]] = []
    if not use_reference_bundle:
        for method in calibration_methods_requested:
            if method == CONTRAST_METHOD_NONE:
                continue
            calibration_methods_skipped.setdefault(
                method,
                "disabled because --use_reference_bundle=false",
            )
        calibration_methods = [m for m in calibration_methods if m == CONTRAST_METHOD_NONE]

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
    )
    need_distal = PROGRAM_DISTAL_ACTIVITY in program_methods or PROGRAM_ENHANCER_BIAS in program_methods
    promoter_links = (
        links_by_method["promoter_overlap"]
        if need_promoter and "promoter_overlap" in links_by_method
        else (_link(peaks, genes, args, "promoter_overlap") if need_promoter else [])
    )
    promoter_indices = promoter_peak_indices(promoter_links) if promoter_links else set()

    needs_ref_ubiquity = (
        PROGRAM_REF_UBIQUITY_PENALTY in calibration_methods
        or CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY in calibration_methods
    )
    needs_atlas = PROGRAM_ATLAS_RESIDUAL in calibration_methods
    manifest_resources: dict[str, dict[str, object]] = {}
    manifest_label = "bundled"
    manifest_warnings: list[str] = []
    ref_peak_idf: list[float] = []
    ref_overlap_count = 0
    atlas_stats_by_definition: dict[str, dict[str, tuple[float, float]]] = {}
    if needs_ref_ubiquity or needs_atlas:
        manifest_label, manifest_resources, _presets, manifest_warnings = load_manifest(_arg(args, "resources_manifest", None))

    if needs_ref_ubiquity:
        ref_resource_id = _arg(args, "ref_ubiquity_resource_id", None) or _default_ref_ubiquity_resource_id(str(args.genome_build))
        resources_dir = (
            Path(_arg(args, "resources_dir", "")).expanduser()
            if _arg(args, "resources_dir", None)
            else default_resources_dir()
        )
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
            calibration_methods_skipped[PROGRAM_REF_UBIQUITY_PENALTY] = str(exc)
            calibration_methods = [m for m in calibration_methods if m != PROGRAM_REF_UBIQUITY_PENALTY]

    if needs_atlas:
        atlas_resource_id = _arg(args, "atlas_resource_id", None) or _default_atlas_resource_id(str(args.genome_build))
        resources_dir = (
            Path(_arg(args, "resources_dir", "")).expanduser()
            if _arg(args, "resources_dir", None)
            else default_resources_dir()
        )
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
            calibration_methods_skipped[PROGRAM_ATLAS_RESIDUAL] = str(exc)
            calibration_methods = [m for m in calibration_methods if m != PROGRAM_ATLAS_RESIDUAL]
    if not calibration_methods:
        calibration_methods = [CONTRAST_METHOD_NONE]
    calibration_methods, auto_reason = resolve_auto_calibration_methods(
        calibration_methods,
        ref_ubiquity_ready=bool(ref_peak_idf),
    )
    if auto_reason:
        calibration_methods_skipped[CONTRAST_METHOD_AUTO_PREFER_REF_UBIQUITY] = auto_reason
    primary_calibration_method = calibration_methods[0]
    _warn_skipped_calibration_methods(
        calibration_methods_skipped,
        str(args.genome_build),
        resource_policy,
    )

    gene_symbol_by_id = {g.gene_id: g.gene_symbol for g in genes}
    gene_biotype_by_id = {g.gene_id: g.gene_biotype for g in genes}

    transforms_by_direction = {"OPEN": "positive", "CLOSE": "negative"}
    full_rows_by_contrast_by_method_by_direction: dict[str, dict[str, dict[str, list[dict[str, object]]]]] = {}
    full_scores_by_contrast_by_method_by_direction: dict[str, dict[str, dict[str, dict[str, float]]]] = {}
    additional_program_rows_by_contrast_by_direction: dict[
        str, dict[str, dict[str, dict[str, list[dict[str, object]]]]]
    ] = {}
    atlas_metric = str(_arg(args, "atlas_metric", "zscore"))
    atlas_eps = float(_arg(args, "atlas_eps", 1e-6))
    atlas_min_raw_quantile = float(_arg(args, "atlas_min_raw_quantile", 0.95))
    atlas_use_log1p = bool(_arg(args, "atlas_use_log1p", True))
    atlas_missing_score_definitions: set[str] = set()
    for direction, transform_mode in transforms_by_direction.items():
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
            return out

        full_scores_none_by_method = _linked_scores_for_peak_values(peak_stat)
        family_scores_none_by_method = _family_scores_for_peak_values(peak_stat)
        full_scores_by_contrast_by_method: dict[str, dict[str, dict[str, float]]] = {}
        family_scores_by_contrast_by_method: dict[str, dict[str, dict[str, dict[str, float]]]] = {}
        for calibration_method in calibration_methods:
            if calibration_method == PROGRAM_ATLAS_RESIDUAL:
                full_scores_by_contrast_by_method[calibration_method] = {}
                family_scores_by_contrast_by_method[calibration_method] = {}
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
                        full_scores_by_contrast_by_method[calibration_method][method] = {}
                    else:
                        atlas_scores = atlas_residual_scores(
                            full_scores_none_by_method[method],
                            linked_stats,
                            atlas_metric,
                            atlas_eps,
                            min_raw_quantile=atlas_min_raw_quantile,
                            use_log1p=atlas_use_log1p,
                        )
                        full_scores_by_contrast_by_method[calibration_method][method] = {
                            g: float(s) for g, s in atlas_scores.items() if float(s) != 0.0
                        }
                    family_scores_by_contrast_by_method[calibration_method][method] = {}
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
                            family_scores_by_contrast_by_method[calibration_method][method][program_method] = {}
                            continue
                        residual_scores = atlas_residual_scores(
                            base_scores,
                            score_stats,
                            atlas_metric,
                            atlas_eps,
                            min_raw_quantile=atlas_min_raw_quantile,
                            use_log1p=atlas_use_log1p,
                        )
                        family_scores_by_contrast_by_method[calibration_method][method][program_method] = {
                            g: float(s) for g, s in residual_scores.items() if float(s) != 0.0
                        }
                continue

            contrast_peak_values = peak_values_for_contrast(
                peak_stat,
                calibration_method,
                ref_peak_idf if ref_peak_idf else None,
            )
            full_scores_by_contrast_by_method[calibration_method] = _linked_scores_for_peak_values(contrast_peak_values)
            family_scores_by_contrast_by_method[calibration_method] = _family_scores_for_peak_values(contrast_peak_values)

        full_scores_by_contrast_by_method_by_direction[direction] = full_scores_by_contrast_by_method

        rows_for_direction: dict[str, dict[str, list[dict[str, object]]]] = {}
        for calibration_method in calibration_methods:
            rows_for_direction[calibration_method] = {}
            for method in link_methods:
                rows_for_direction[calibration_method][method] = _rows_from_scores(
                    full_scores_by_contrast_by_method[calibration_method][method],
                    gene_symbol_by_id,
                    gene_biotype_by_id,
                )
        full_rows_by_contrast_by_method_by_direction[direction] = rows_for_direction

        program_rows_by_contrast_by_method: dict[str, dict[str, dict[str, list[dict[str, object]]]]] = {}
        for calibration_method in calibration_methods:
            rows_by_method: dict[str, dict[str, list[dict[str, object]]]] = {}
            for method in link_methods:
                rows_by_program: dict[str, list[dict[str, object]]] = {}
                for program_method in program_methods:
                    if program_method == PROGRAM_LINKED_ACTIVITY:
                        continue
                    rows_by_program[program_method] = _rows_from_scores(
                        family_scores_by_contrast_by_method.get(calibration_method, {})
                        .get(method, {})
                        .get(program_method, {}),
                        gene_symbol_by_id,
                        gene_biotype_by_id,
                    )
                rows_by_method[method] = rows_by_program
            program_rows_by_contrast_by_method[calibration_method] = rows_by_method
        additional_program_rows_by_contrast_by_direction[direction] = program_rows_by_contrast_by_method

    if atlas_missing_score_definitions:
        msg = (
            "atlas_residual baseline missing for score definitions: "
            + ", ".join(sorted(atlas_missing_score_definitions))
        )
        calibration_methods_skipped.setdefault(PROGRAM_ATLAS_RESIDUAL, msg)
        print(f"warning: {msg}", file=sys.stderr)
        print(
            "warning: provide per-score-definition atlas baselines or run with --calibration_methods none/ref_ubiquity_penalty.",
            file=sys.stderr,
        )

    selected_scores = full_scores_by_contrast_by_method_by_direction[selected_direction][primary_calibration_method][primary_link_method]
    selected_gene_ids = _select_gene_ids(selected_scores, args)
    selected_weights = _selected_weights(selected_scores, selected_gene_ids, normalize)

    selected_rows: list[dict[str, object]] = []
    for rank, gene_id in enumerate(selected_gene_ids, start=1):
        selected_rows.append(
            {
                "gene_id": gene_id,
                "score": float(selected_scores[gene_id]),
                "weight": float(selected_weights.get(gene_id, 0.0)),
                "rank": rank,
                "gene_symbol": gene_symbol_by_id.get(gene_id),
                "gene_biotype": gene_biotype_by_id.get(gene_id),
            }
        )
    _write_rows(out_dir / "geneset.tsv", selected_rows)

    full_rows = full_rows_by_contrast_by_method_by_direction[selected_direction][primary_calibration_method][primary_link_method]
    if emit_full:
        _write_rows(out_dir / "geneset.full.tsv", full_rows)

    gmt_sets: list[tuple[str, list[str]]] = []
    gmt_plans: list[dict[str, object]] = []
    requested_gmt_outputs: list[dict[str, str]] = []
    gmt_diagnostics: list[dict[str, object]] = []
    gmt_path = resolve_gmt_out_path(out_dir, gmt_out)
    linked_output_methods = link_methods_for_program(program_preset, PROGRAM_LINKED_ACTIVITY, link_methods)
    if PROGRAM_LINKED_ACTIVITY in program_methods and not linked_output_methods:
        program_methods_skipped[PROGRAM_LINKED_ACTIVITY] = (
            "no compatible link_method selected for preset; expected nearest_tss for connectable/default"
        )
    if emit_gmt:
        for direction in ("OPEN", "CLOSE"):
            base_prefix = (
                f"atac_bulk_matrix__dataset={dataset_label}"
                "__study_contrast=condition_between_samples"
                f"__condition_column={condition_column}"
                f"__condition_a={condition_a}"
                f"__condition_b={condition_b}"
                f"__direction={direction}"
            )
            if PROGRAM_LINKED_ACTIVITY in program_methods:
                for calibration_method in calibration_methods:
                    for method in linked_output_methods:
                        requested_gmt_outputs.append(
                            {
                                "program_method": PROGRAM_LINKED_ACTIVITY,
                                "calibration_method": calibration_method,
                                "link_method": method,
                                "direction": direction,
                            }
                        )
                        method_sets, method_plans = build_gmt_sets_from_rows(
                            rows=full_rows_by_contrast_by_method_by_direction[direction][calibration_method][method],
                            base_name=(
                                f"{base_prefix}__program={PROGRAM_LINKED_ACTIVITY}"
                                f"__calibration_method={calibration_method}"
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
                            diagnostics=gmt_diagnostics,
                            context={
                                "program_method": PROGRAM_LINKED_ACTIVITY,
                                "calibration_method": calibration_method,
                                "link_method": method,
                                "direction": direction,
                            },
                        )
                        for plan in method_plans:
                            params = dict(plan.get("parameters", {}))
                            params["program_method"] = PROGRAM_LINKED_ACTIVITY
                            params["calibration_method"] = calibration_method
                            params["link_method"] = method
                            params["direction"] = direction
                            plan["parameters"] = params
                        gmt_sets.extend(method_sets)
                        gmt_plans.extend(method_plans)

            for calibration_method, rows_by_method in additional_program_rows_by_contrast_by_direction[direction].items():
                for method, rows_by_program in rows_by_method.items():
                    for program_method, rows in rows_by_program.items():
                        allowed_methods = link_methods_for_program(program_preset, program_method, link_methods)
                        if method not in allowed_methods:
                            continue
                        requested_gmt_outputs.append(
                            {
                                "program_method": program_method,
                                "calibration_method": calibration_method,
                                "link_method": method,
                                "direction": direction,
                            }
                        )
                        method_sets, method_plans = build_gmt_sets_from_rows(
                            rows=rows,
                            base_name=(
                                f"{base_prefix}__program={program_method}"
                                f"__calibration_method={calibration_method}"
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
                            diagnostics=gmt_diagnostics,
                            context={
                                "program_method": program_method,
                                "calibration_method": calibration_method,
                                "link_method": method,
                                "direction": direction,
                            },
                        )
                        for plan in method_plans:
                            params = dict(plan.get("parameters", {}))
                            params["program_method"] = program_method
                            params["calibration_method"] = calibration_method
                            params["link_method"] = method
                            params["direction"] = direction
                            plan["parameters"] = params
                        gmt_sets.extend(method_sets)
                        gmt_plans.extend(method_plans)
        for program_method in program_methods:
            if program_method == PROGRAM_LINKED_ACTIVITY:
                continue
            allowed_methods = link_methods_for_program(program_preset, program_method, link_methods)
            if allowed_methods:
                continue
            if program_method not in program_methods_skipped:
                program_methods_skipped[program_method] = "no compatible link_method selected for active preset"
        write_gmt(gmt_sets, gmt_path)
    _warn_gmt_output_diagnostics(gmt_diagnostics)
    _warn_skipped_program_methods(program_methods_skipped)

    resources_dir = (
        Path(_arg(args, "resources_dir", "")).expanduser()
        if _arg(args, "resources_dir", None)
        else default_resources_dir()
    )
    files = [
        input_file_record(args.peak_bed, "peak_bed"),
        input_file_record(args.peak_matrix_tsv, "peak_matrix_tsv"),
        input_file_record(args.sample_metadata_tsv, "sample_metadata_tsv"),
        input_file_record(args.gtf, "gtf"),
    ]
    if _arg(args, "qc_marker_genes_tsv", None):
        files.append(input_file_record(_arg(args, "qc_marker_genes_tsv", None), "qc_marker_genes_tsv"))
    if _EXTERNAL_LINK_METHOD in link_methods:
        files.append(input_file_record(args.region_gene_links_tsv, "region_gene_links_tsv"))
    for r in resources_used:
        files.append(input_file_record(str(r["path"]), f"resource:{r['id']}"))

    assigned_primary = len({int(link["peak_index"]) for link in links_by_method[primary_link_method]})
    link_assignment: dict[str, dict[str, float]] = {}
    for method in link_methods:
        assigned = len({int(link["peak_index"]) for link in links_by_method[method]})
        link_assignment[method] = {
            "n_features_assigned": int(assigned),
            "fraction_features_assigned": float(assigned / len(peaks) if peaks else 0.0),
        }
    marker_qc = marker_hit_summary(selected_rows, marker_genes)
    emitted_combinations = collect_emitted_method_combinations(gmt_plans)
    skipped_gmt_outputs = _collect_skipped_gmt_outputs(gmt_diagnostics)
    if not emitted_combinations:
        emitted_combinations = [
            {
                "program_method": PROGRAM_LINKED_ACTIVITY,
                "calibration_method": primary_calibration_method,
                "link_method": primary_link_method,
                "direction": selected_direction,
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

    output_files = [{"path": str(out_dir / "geneset.tsv"), "role": "selected_program"}]
    if emit_full:
        output_files.append({"path": str(out_dir / "geneset.full.tsv"), "role": "full_scores"})
    if emit_gmt:
        output_files.append({"path": str(gmt_path), "role": "gmt"})
    run_summary_payload: dict[str, object] = {
        "converter": "atac_bulk_matrix",
        "dataset_label": dataset_label,
        "program_preset": program_preset,
        "primary_program_method": PROGRAM_LINKED_ACTIVITY,
        "primary_link_method": primary_link_method,
        "study_contrast": "condition_between_samples",
        "primary_calibration_method": primary_calibration_method,
        "selected_direction": selected_direction,
        "n_input_peaks": len(peaks),
        "link_assignment": link_assignment,
        "promoter_peak_count": len(promoter_indices),
        "distal_peak_count": len(peaks) - len(promoter_indices),
        "emitted_method_combinations": emitted_combinations,
        "requested_gmt_outputs": requested_gmt_outputs,
        "skipped_gmt_outputs": skipped_gmt_outputs,
        "program_methods_skipped": program_methods_skipped,
        "calibration_methods_skipped": calibration_methods_skipped,
    }
    if marker_qc is not None:
        run_summary_payload["marker_qc"] = marker_qc
    if ref_summary is not None:
        run_summary_payload["ref_ubiquity"] = ref_summary
    run_summary_json_path, run_summary_txt_path = write_run_summary_files(out_dir, run_summary_payload)
    output_files.append({"path": str(run_summary_json_path), "role": "run_summary_json"})
    output_files.append({"path": str(run_summary_txt_path), "role": "run_summary_text"})

    params = _resolved_parameters(args, condition_a, condition_b)
    params["program_methods_active"] = program_methods
    params["program_methods_skipped"] = program_methods_skipped
    params["calibration_methods_active"] = calibration_methods
    params["calibration_methods_skipped"] = calibration_methods_skipped
    params["primary_calibration_method"] = primary_calibration_method
    params["selected_direction"] = selected_direction
    summary_payload: dict[str, object] = {
        "n_input_features": len(peaks),
        "n_genes": len(selected_rows),
        "n_features_assigned": assigned_primary,
        "fraction_features_assigned": (assigned_primary / len(peaks) if peaks else 0.0),
        "n_external_links_unresolved_gene_id": external_links_unresolved,
        "n_program_methods": len(program_methods),
        "n_calibration_methods": len(calibration_methods),
        "n_calibration_methods_skipped": len(calibration_methods_skipped),
        "n_resource_manifest_warnings": len(manifest_warnings),
        "n_samples": len(sample_ids),
        "n_samples_a": len(a_indices),
        "n_samples_b": len(b_indices),
    }
    if marker_qc is not None:
        summary_payload["marker_qc"] = marker_qc
    meta = make_metadata(
        converter_name="atac_bulk_matrix",
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
            "weight_type": "signed",
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
            "score_definition": "sum over peaks of transformed differential peak_stat times link_weight",
            "study_contrast": {
                "mode": "condition_between_samples",
                "contrast_type": "condition_between_samples",
                "condition_column": condition_column,
                "condition_a": condition_a,
                "condition_b": condition_b,
                "metric": contrast_metric,
                "pseudocount": pseudocount,
                "n_samples": len(sample_ids),
                "n_samples_a": len(a_indices),
                "n_samples_b": len(b_indices),
                "directions_emitted": ["OPEN", "CLOSE"],
                "selected_direction": selected_direction,
            },
            "contrast": {
                "mode": "condition_between_samples",
                "contrast_type": "condition_between_samples",
                "condition_column": condition_column,
                "condition_a": condition_a,
                "condition_b": condition_b,
                "metric": contrast_metric,
                "pseudocount": pseudocount,
                "n_samples": len(sample_ids),
                "n_samples_a": len(a_indices),
                "n_samples_b": len(b_indices),
                "directions_emitted": ["OPEN", "CLOSE"],
                "selected_direction": selected_direction,
            },
            "calibration_methods": calibration_methods,
            "calibration_methods_skipped": calibration_methods_skipped,
            "primary_calibration_method": primary_calibration_method,
            "program_methods": program_methods,
            "program_methods_skipped": program_methods_skipped,
            "primary_program_method": PROGRAM_LINKED_ACTIVITY,
        },
        output_files=output_files,
        gmt={
            "written": emit_gmt,
            "path": (str(gmt_path.relative_to(out_dir)) if gmt_path.is_relative_to(out_dir) else str(gmt_path))
            if emit_gmt
            else None,
            "prefer_symbol": gmt_prefer_symbol,
            "require_symbol": gmt_require_symbol,
            "biotype_allowlist": gmt_biotype_allowlist,
            "min_genes": gmt_min_genes,
            "max_genes": gmt_max_genes,
            "emit_small_gene_sets": emit_small_gene_sets,
            "requested_outputs": requested_gmt_outputs,
            "emitted_outputs": emitted_combinations,
            "skipped_outputs": skipped_gmt_outputs,
            "diagnostics": gmt_diagnostics,
            "plans": gmt_plans,
        },
    )
    if resources_used or manifest_warnings:
        meta["resources"] = {
            "manifest": manifest_label,
            "resources_dir": str(resources_dir),
            "used": resources_used,
            "warnings": manifest_warnings,
        }
    write_metadata(out_dir / "geneset.meta.json", meta)

    return {
        "n_peaks": len(peaks),
        "n_genes": len(selected_rows),
        "out_dir": str(out_dir),
        "program_methods": program_methods,
        "program_methods_skipped": program_methods_skipped,
        "calibration_methods": calibration_methods,
        "calibration_methods_skipped": calibration_methods_skipped,
    }
