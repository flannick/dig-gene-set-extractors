from __future__ import annotations

import csv
import math
from pathlib import Path

from omics2geneset.core.atac_programs import (
    PROGRAM_ATLAS_RESIDUAL,
    PROGRAM_DISTAL_ACTIVITY,
    PROGRAM_ENHANCER_BIAS,
    PROGRAM_LINKED_ACTIVITY,
    PROGRAM_PROMOTER_ACTIVITY,
    PROGRAM_REF_UBIQUITY_PENALTY,
    atlas_residual_scores,
    enhancer_bias_scores,
    mask_peak_weights,
    promoter_peak_indices,
    resolve_program_methods,
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
from omics2geneset.core.peak_to_gene import (
    link_distance_decay,
    link_external_regions,
    link_nearest_tss,
    link_promoter_overlap,
)
from omics2geneset.core.reference_calibration import apply_peak_idf, peak_ref_idf_by_overlap
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
from omics2geneset.io.reference_tables import read_atlas_gene_stats_tsv, read_ref_ubiquity_tsv
from omics2geneset.io.region_gene_links import read_region_gene_links_tsv
from omics2geneset.resource_manager import default_resources_dir, load_manifest, resource_metadata_record


_LINK_METHODS = ("promoter_overlap", "nearest_tss", "distance_decay")
_EXTERNAL_LINK_METHOD = "external"


def _arg(args, name: str, default):
    return getattr(args, name, default)


def _resolve_dataset_label(args) -> str:
    label = _arg(args, "dataset_label", None)
    if label is None:
        return "dataset"
    text = str(label).strip()
    return text or "dataset"


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


def _default_ref_ubiquity_resource_id(genome_build: str) -> str | None:
    gb = genome_build.strip().lower()
    if gb in {"hg38", "grch38"}:
        return "ccre_ubiquity_hg38"
    if gb in {"mm10", "grcm38"}:
        return "ccre_ubiquity_mm10"
    return None


def _default_atlas_resource_id(genome_build: str) -> str | None:
    gb = genome_build.strip().lower()
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
    gmt_topk_list = parse_int_list_csv(str(_arg(args, "gmt_topk_list", "100,200,500")))
    gmt_mass_list = parse_mass_list_csv(str(_arg(args, "gmt_mass_list", "0.5,0.8,0.9")))
    gmt_biotype_allowlist = parse_str_list_csv(str(_arg(args, "gmt_biotype_allowlist", "protein_coding")))
    link_methods = _resolve_link_methods(args)
    program_methods = resolve_program_methods(
        "bulk",
        _arg(args, "program_preset", "connectable"),
        _arg(args, "program_methods", None),
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
        "program_preset": _arg(args, "program_preset", "connectable"),
        "program_methods": _arg(args, "program_methods", None),
        "program_methods_evaluated": program_methods,
        "resource_policy": _arg(args, "resource_policy", "skip"),
        "ref_ubiquity_resource_id": ref_ubiquity_resource_id,
        "atlas_resource_id": atlas_resource_id,
        "atlas_metric": _arg(args, "atlas_metric", "logratio"),
        "atlas_eps": _arg(args, "atlas_eps", 1e-6),
        "promoter_upstream_bp": _arg(args, "promoter_upstream_bp", 2000),
        "promoter_downstream_bp": _arg(args, "promoter_downstream_bp", 500),
        "max_distance_bp": max_distance,
        "decay_length_bp": _arg(args, "decay_length_bp", 50000),
        "max_genes_per_peak": _arg(args, "max_genes_per_peak", 5),
        "external_linking_enabled": _EXTERNAL_LINK_METHOD in link_methods,
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
    dataset_label = _resolve_dataset_label(args)

    resource_policy = str(_arg(args, "resource_policy", "skip"))
    if resource_policy not in {"skip", "fail"}:
        raise ValueError(f"Unsupported resource_policy: {resource_policy}")
    program_methods = resolve_program_methods(
        "bulk",
        _arg(args, "program_preset", "connectable"),
        _arg(args, "program_methods", None),
    )
    program_methods_skipped: dict[str, str] = {}
    resources_used: list[dict[str, object]] = []

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
    distance_links = (
        links_by_method["distance_decay"]
        if need_distal and "distance_decay" in links_by_method
        else (_link(peaks, genes, args, "distance_decay") if need_distal else [])
    )
    promoter_indices = promoter_peak_indices(promoter_links) if promoter_links else set()

    needs_ref_ubiquity = PROGRAM_REF_UBIQUITY_PENALTY in program_methods
    needs_atlas = PROGRAM_ATLAS_RESIDUAL in program_methods
    manifest_resources: dict[str, dict[str, object]] = {}
    manifest_label = "bundled"
    manifest_warnings: list[str] = []
    ref_peak_idf: list[float] = []
    atlas_stats: dict[str, tuple[float, float]] = {}
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
            program_methods_skipped[PROGRAM_REF_UBIQUITY_PENALTY] = str(exc)
            program_methods = [m for m in program_methods if m != PROGRAM_REF_UBIQUITY_PENALTY]

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
            atlas_stats = read_atlas_gene_stats_tsv(atlas_path)
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
            program_methods_skipped[PROGRAM_ATLAS_RESIDUAL] = str(exc)
            program_methods = [m for m in program_methods if m != PROGRAM_ATLAS_RESIDUAL]

    gene_symbol_by_id = {g.gene_id: g.gene_symbol for g in genes}
    gene_biotype_by_id = {g.gene_id: g.gene_biotype for g in genes}

    transforms_by_direction = {"OPEN": "positive", "CLOSE": "negative"}
    full_rows_by_method_by_direction: dict[str, dict[str, list[dict[str, object]]]] = {}
    full_scores_by_method_by_direction: dict[str, dict[str, dict[str, float]]] = {}
    additional_program_rows_by_direction: dict[str, dict[str, list[dict[str, object]]]] = {}
    for direction, transform_mode in transforms_by_direction.items():
        full_scores_by_method: dict[str, dict[str, float]] = {}
        for method in link_methods:
            raw_scores = score_genes(peak_stat, links_by_method[method], transform_mode)
            full_scores_by_method[method] = {g: float(s) for g, s in raw_scores.items() if float(s) != 0.0}
        full_scores_by_method_by_direction[direction] = full_scores_by_method

        rows_for_direction: dict[str, list[dict[str, object]]] = {}
        for method in link_methods:
            rows_for_direction[method] = _rows_from_scores(
                full_scores_by_method[method],
                gene_symbol_by_id,
                gene_biotype_by_id,
            )
        full_rows_by_method_by_direction[direction] = rows_for_direction

        program_rows: dict[str, list[dict[str, object]]] = {}
        promoter_scores: dict[str, float] = {}
        distal_scores: dict[str, float] = {}
        if PROGRAM_PROMOTER_ACTIVITY in program_methods or PROGRAM_ENHANCER_BIAS in program_methods:
            promoter_peak_stat = mask_peak_weights(peak_stat, include_indices=promoter_indices)
            promoter_raw = score_genes(promoter_peak_stat, promoter_links, transform_mode)
            promoter_scores = {g: float(s) for g, s in promoter_raw.items() if float(s) != 0.0}
            if PROGRAM_PROMOTER_ACTIVITY in program_methods:
                program_rows[PROGRAM_PROMOTER_ACTIVITY] = _rows_from_scores(
                    promoter_scores,
                    gene_symbol_by_id,
                    gene_biotype_by_id,
                )
        if PROGRAM_DISTAL_ACTIVITY in program_methods or PROGRAM_ENHANCER_BIAS in program_methods:
            distal_peak_stat = mask_peak_weights(peak_stat, exclude_indices=promoter_indices)
            distal_raw = score_genes(distal_peak_stat, distance_links, transform_mode)
            distal_scores = {g: float(s) for g, s in distal_raw.items() if float(s) != 0.0}
            if PROGRAM_DISTAL_ACTIVITY in program_methods:
                program_rows[PROGRAM_DISTAL_ACTIVITY] = _rows_from_scores(
                    distal_scores,
                    gene_symbol_by_id,
                    gene_biotype_by_id,
                )
        if PROGRAM_ENHANCER_BIAS in program_methods:
            bias_scores = enhancer_bias_scores(promoter_scores, distal_scores)
            bias_nonzero = {g: float(s) for g, s in bias_scores.items() if float(s) != 0.0}
            program_rows[PROGRAM_ENHANCER_BIAS] = _rows_from_scores(
                bias_nonzero,
                gene_symbol_by_id,
                gene_biotype_by_id,
            )
        if PROGRAM_REF_UBIQUITY_PENALTY in program_methods:
            ref_peak_stat = apply_peak_idf(peak_stat, ref_peak_idf)
            ref_raw = score_genes(ref_peak_stat, links_by_method[primary_link_method], transform_mode)
            ref_scores = {g: float(s) for g, s in ref_raw.items() if float(s) != 0.0}
            program_rows[PROGRAM_REF_UBIQUITY_PENALTY] = _rows_from_scores(
                ref_scores,
                gene_symbol_by_id,
                gene_biotype_by_id,
            )
        if PROGRAM_ATLAS_RESIDUAL in program_methods:
            atlas_metric = str(_arg(args, "atlas_metric", "logratio"))
            atlas_eps = float(_arg(args, "atlas_eps", 1e-6))
            atlas_scores = atlas_residual_scores(
                full_scores_by_method[primary_link_method],
                atlas_stats,
                atlas_metric,
                atlas_eps,
            )
            atlas_nonzero = {g: float(s) for g, s in atlas_scores.items() if float(s) != 0.0}
            program_rows[PROGRAM_ATLAS_RESIDUAL] = _rows_from_scores(
                atlas_nonzero,
                gene_symbol_by_id,
                gene_biotype_by_id,
            )
        additional_program_rows_by_direction[direction] = program_rows

    selected_scores = score_genes(peak_stat, links_by_method[primary_link_method], selected_transform)
    selected_scores = {g: float(s) for g, s in selected_scores.items() if float(s) != 0.0}
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

    full_rows = full_rows_by_method_by_direction[selected_direction][primary_link_method]
    if emit_full:
        _write_rows(out_dir / "geneset.full.tsv", full_rows)

    gmt_sets: list[tuple[str, list[str]]] = []
    gmt_plans: list[dict[str, object]] = []
    gmt_path = resolve_gmt_out_path(out_dir, gmt_out)
    if emit_gmt:
        for direction in ("OPEN", "CLOSE"):
            base_prefix = (
                f"atac_bulk_matrix__dataset={dataset_label}"
                "__contrast=condition_between_samples"
                f"__condition_column={condition_column}"
                f"__condition_a={condition_a}"
                f"__condition_b={condition_b}"
                f"__direction={direction}"
            )
            if PROGRAM_LINKED_ACTIVITY in program_methods:
                for method in link_methods:
                    method_sets, method_plans = build_gmt_sets_from_rows(
                        rows=full_rows_by_method_by_direction[direction][method],
                        base_name=f"{base_prefix}__program={PROGRAM_LINKED_ACTIVITY}__link_method={method}",
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
                        params["direction"] = direction
                        plan["parameters"] = params
                    gmt_sets.extend(method_sets)
                    gmt_plans.extend(method_plans)

            for program_method, rows in additional_program_rows_by_direction[direction].items():
                method_sets, method_plans = build_gmt_sets_from_rows(
                    rows=rows,
                    base_name=f"{base_prefix}__program={program_method}",
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
                    params["direction"] = direction
                    plan["parameters"] = params
                gmt_sets.extend(method_sets)
                gmt_plans.extend(method_plans)
        write_gmt(gmt_sets, gmt_path)

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
    if _EXTERNAL_LINK_METHOD in link_methods:
        files.append(input_file_record(args.region_gene_links_tsv, "region_gene_links_tsv"))
    for r in resources_used:
        files.append(input_file_record(str(r["path"]), f"resource:{r['id']}"))

    output_files = [{"path": str(out_dir / "geneset.tsv"), "role": "selected_program"}]
    if emit_full:
        output_files.append({"path": str(out_dir / "geneset.full.tsv"), "role": "full_scores"})
    if emit_gmt:
        output_files.append({"path": str(gmt_path), "role": "gmt"})

    params = _resolved_parameters(args, condition_a, condition_b)
    params["program_methods_active"] = program_methods
    params["program_methods_skipped"] = program_methods_skipped
    params["selected_direction"] = selected_direction
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
        summary={
            "n_input_features": len(peaks),
            "n_genes": len(selected_rows),
            "n_features_assigned": len({int(link["peak_index"]) for link in links_by_method[primary_link_method]}),
            "fraction_features_assigned": (
                len({int(link["peak_index"]) for link in links_by_method[primary_link_method]}) / len(peaks)
                if peaks
                else 0.0
            ),
            "n_external_links_unresolved_gene_id": external_links_unresolved,
            "n_program_methods": len(program_methods),
            "n_resource_methods_skipped": len(program_methods_skipped),
            "n_resource_manifest_warnings": len(manifest_warnings),
            "n_samples": len(sample_ids),
            "n_samples_a": len(a_indices),
            "n_samples_b": len(b_indices),
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
            "score_definition": "sum over peaks of transformed differential peak_stat times link_weight",
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
    }
