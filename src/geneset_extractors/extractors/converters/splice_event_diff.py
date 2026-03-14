from __future__ import annotations

import csv
import os
from pathlib import Path

from geneset_extractors.core.metadata import input_file_record
from geneset_extractors.extractors.splicing.resource_utils import (
    build_resources_info,
    load_resource_context,
    resolve_resource_path,
)
from geneset_extractors.extractors.splicing.splice_event_diff_workflow import (
    SpliceRow,
    SplicingWorkflowConfig,
    read_event_alias_tsv,
    read_event_impact_tsv,
    read_event_ubiquity_tsv,
    read_tsv_rows,
    run_splice_event_diff_workflow,
)


DEFAULT_ALIAS_RESOURCE_ID = "splice_event_aliases_human_v1"
DEFAULT_UBIQUITY_RESOURCE_ID = "splice_event_ubiquity_human_v1"
DEFAULT_IMPACT_RESOURCE_ID = "splice_event_impact_human_v1"



def _default_resource_id(organism: str, kind: str) -> str | None:
    if str(organism).strip().lower() != "human":
        return None
    if kind == "alias":
        return DEFAULT_ALIAS_RESOURCE_ID
    if kind == "ubiquity":
        return DEFAULT_UBIQUITY_RESOURCE_ID
    if kind == "impact":
        return DEFAULT_IMPACT_RESOURCE_ID
    return None



def _clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()



def _merge_leafcutter_cluster_stats(rows: list[SpliceRow], cluster_stats_tsv: str, event_group_column: str | None) -> tuple[list[SpliceRow], bool]:
    if not cluster_stats_tsv:
        return rows, False
    cluster_col = event_group_column or "cluster"
    with Path(cluster_stats_tsv).open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"LeafCutter cluster stats table has no header: {cluster_stats_tsv}")
        candidate_cluster_cols = [cluster_col, "cluster", "cluster_id"]
        candidate_p_cols = ["padj", "fdr", "p.adjust", "pvalue", "pval"]
        resolved_cluster = next((c for c in candidate_cluster_cols if c in reader.fieldnames), None)
        resolved_p = next((c for c in candidate_p_cols if c in reader.fieldnames), None)
        if resolved_cluster is None or resolved_p is None:
            return rows, False
        cluster_map: dict[str, str] = {}
        for row in reader:
            key = _clean(row.get(resolved_cluster))
            value = _clean(row.get(resolved_p))
            if key and value:
                cluster_map[key] = value
    if not cluster_map:
        return rows, False

    merged: list[SpliceRow] = []
    used = False
    for row in rows:
        values = dict(row.values)
        cluster_id = _clean(values.get(cluster_col)) or _clean(values.get("cluster")) or _clean(values.get("cluster_id"))
        if cluster_id and cluster_id in cluster_map:
            values.setdefault("padj", cluster_map[cluster_id])
            values.setdefault("pvalue", cluster_map[cluster_id])
            used = True
        merged.append(SpliceRow(line_no=row.line_no, values=values))
    return merged, used



def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fieldnames, rows = read_tsv_rows(args.splice_tsv)
    resource_policy = str(args.resource_policy)
    if resource_policy not in {"skip", "fail"}:
        raise ValueError(f"Unsupported resource_policy: {resource_policy}")

    has_resource_location = bool(
        args.resources_dir
        or args.resources_manifest
        or os.getenv("GENESET_EXTRACTORS_RESOURCES_DIR")
        or os.getenv("OMICS2GENESET_RESOURCES_DIR")
    )
    use_bundle = bool(args.use_reference_bundle)
    need_resources = bool(
        use_bundle
        and (
            has_resource_location
            or args.event_alias_resource_id
            or args.event_ubiquity_resource_id
            or args.event_impact_resource_id
        )
    )
    ctx = load_resource_context(args.resources_manifest, args.resources_dir) if need_resources else None

    alias_path = None
    alias_resource_id = args.event_alias_resource_id or _default_resource_id(args.organism, "alias")
    if use_bundle and alias_resource_id and ctx is not None:
        alias_path = resolve_resource_path(
            ctx=ctx,
            resource_id=alias_resource_id,
            resource_policy=resource_policy,
            role_label="event_alias",
            enablement_hint=(
                "Provide --resources_dir containing the splice-event alias table, or set --event_alias_resource_id "
                "to a valid local resource."
            ),
        )

    ubiquity_path = None
    ubiquity_resource_id = args.event_ubiquity_resource_id or _default_resource_id(args.organism, "ubiquity")
    if use_bundle and ubiquity_resource_id and ctx is not None:
        ubiquity_path = resolve_resource_path(
            ctx=ctx,
            resource_id=ubiquity_resource_id,
            resource_policy=resource_policy,
            role_label="event_ubiquity",
            enablement_hint=(
                "Provide --resources_dir containing the splice-event ubiquity table, or set --event_ubiquity_resource_id "
                "to a valid local resource."
            ),
        )

    impact_path = None
    impact_resource_id = args.event_impact_resource_id or _default_resource_id(args.organism, "impact")
    if use_bundle and impact_resource_id and ctx is not None:
        impact_path = resolve_resource_path(
            ctx=ctx,
            resource_id=impact_resource_id,
            resource_policy=resource_policy,
            role_label="event_impact",
            enablement_hint=(
                "Provide --resources_dir containing the splice-event impact table, or set --event_impact_resource_id "
                "to a valid local resource."
            ),
        )

    if args.tool_family == "leafcutter" and args.cluster_stats_tsv:
        rows, cluster_stats_used = _merge_leafcutter_cluster_stats(rows, args.cluster_stats_tsv, args.event_group_column)
    else:
        cluster_stats_used = False

    alias_map = read_event_alias_tsv(alias_path) if alias_path is not None else None
    ubiquity_map = read_event_ubiquity_tsv(ubiquity_path) if ubiquity_path is not None else None
    impact_map = read_event_impact_tsv(impact_path) if impact_path is not None else None

    files = [input_file_record(args.splice_tsv, "splice_tsv")]
    if args.cluster_stats_tsv:
        files.append(input_file_record(args.cluster_stats_tsv, "cluster_stats_tsv"))
    if alias_path is not None:
        files.append(input_file_record(alias_path, "event_alias_table"))
    if ubiquity_path is not None:
        files.append(input_file_record(ubiquity_path, "event_ubiquity_table"))
    if impact_path is not None:
        files.append(input_file_record(impact_path, "event_impact_table"))

    dataset_label = str(args.dataset_label or "").strip() or Path(args.splice_tsv).name
    signature_name = str(args.signature_name or "").strip() or Path(args.splice_tsv).stem
    cfg = SplicingWorkflowConfig(
        converter_name="splice_event_diff",
        out_dir=out_dir,
        organism=args.organism,
        genome_build=args.genome_build,
        dataset_label=dataset_label,
        signature_name=signature_name,
        tool_family=args.tool_family,
        event_id_column=args.event_id_column,
        event_group_column=args.event_group_column,
        event_type_column=args.event_type_column,
        gene_id_column=args.gene_id_column,
        gene_symbol_column=args.gene_symbol_column,
        chrom_column=args.chrom_column,
        start_column=args.start_column,
        end_column=args.end_column,
        strand_column=args.strand_column,
        score_column=args.score_column,
        stat_column=args.stat_column,
        delta_psi_column=args.delta_psi_column,
        psi_column=args.psi_column,
        padj_column=args.padj_column,
        pvalue_column=args.pvalue_column,
        probability_column=args.probability_column,
        read_support_column=args.read_support_column,
        novel_flag_column=args.novel_flag_column,
        annotation_status_column=args.annotation_status_column,
        score_mode=args.score_mode,
        score_transform=args.score_transform,
        confidence_weight_mode=args.confidence_weight_mode,
        min_probability=float(args.min_probability),
        min_read_support=float(args.min_read_support),
        neglog10p_cap=float(args.neglog10p_cap),
        neglog10p_eps=float(args.neglog10p_eps),
        event_dup_policy=args.event_dup_policy,
        gene_aggregation=args.gene_aggregation,
        gene_topk_events=int(args.gene_topk_events),
        ambiguous_gene_policy=args.ambiguous_gene_policy,
        impact_mode=args.impact_mode,
        impact_min=float(args.impact_min),
        impact_max=float(args.impact_max),
        select=args.select,
        top_k=int(args.top_k),
        quantile=float(args.quantile),
        min_score=float(args.min_score),
        normalize=args.normalize,
        emit_full=bool(args.emit_full),
        emit_gmt=bool(args.emit_gmt),
        gmt_out=args.gmt_out,
        gmt_format=args.gmt_format,
        gmt_prefer_symbol=bool(args.gmt_prefer_symbol),
        gmt_require_symbol=bool(args.gmt_require_symbol),
        gmt_biotype_allowlist=args.gmt_biotype_allowlist,
        gmt_min_genes=int(args.gmt_min_genes),
        gmt_max_genes=int(args.gmt_max_genes),
        gmt_topk_list=args.gmt_topk_list,
        gmt_mass_list=args.gmt_mass_list,
        gmt_split_signed=bool(args.gmt_split_signed),
        emit_small_gene_sets=bool(args.emit_small_gene_sets),
    )
    resources_info = build_resources_info(ctx) if ctx is not None else None
    result = run_splice_event_diff_workflow(
        cfg=cfg,
        fieldnames=fieldnames,
        rows=rows,
        alias_map=alias_map,
        ubiquity_map=ubiquity_map,
        impact_map=impact_map,
        input_files=files,
        resources_info=resources_info,
    )
    result["cluster_stats_used"] = cluster_stats_used
    return result
