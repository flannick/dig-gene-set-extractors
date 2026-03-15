from __future__ import annotations

import os
from pathlib import Path

from geneset_extractors.core.metadata import input_file_record
from geneset_extractors.extractors.splicing.resource_utils import (
    build_resources_info,
    load_resource_context,
    resolve_resource_path,
)
from geneset_extractors.extractors.splicing.splice_event_diff_workflow import (
    read_event_alias_tsv,
    read_gene_burden_tsv,
    read_event_impact_tsv,
    read_event_ubiquity_tsv,
)
from geneset_extractors.extractors.splicing.splice_event_matrix_workflow import (
    SpliceMatrixWorkflowConfig,
    run_splice_event_matrix_workflow,
)


DEFAULT_ALIAS_RESOURCE_ID = "splice_event_aliases_human_v1"
DEFAULT_UBIQUITY_RESOURCE_ID = "splice_event_ubiquity_human_v1"
DEFAULT_IMPACT_RESOURCE_ID = "splice_event_impact_human_v1"
DEFAULT_GENE_BURDEN_RESOURCE_ID = "splice_gene_event_burden_human_v1"



def _default_resource_id(organism: str, kind: str) -> str | None:
    if str(organism).strip().lower() != "human":
        return None
    if kind == "alias":
        return DEFAULT_ALIAS_RESOURCE_ID
    if kind == "ubiquity":
        return DEFAULT_UBIQUITY_RESOURCE_ID
    if kind == "impact":
        return DEFAULT_IMPACT_RESOURCE_ID
    if kind == "gene_burden":
        return DEFAULT_GENE_BURDEN_RESOURCE_ID
    return None



def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

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
            or getattr(args, "gene_burden_resource_id", None)
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
                "Provide --resources_dir containing the splice-event alias table, or set --event_alias_resource_id to a valid local resource."
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
                "Provide --resources_dir containing the splice-event ubiquity table, or set --event_ubiquity_resource_id to a valid local resource."
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
                "Provide --resources_dir containing the splice-event impact table, or set --event_impact_resource_id to a valid local resource."
            ),
        )
    gene_burden_path = None
    gene_burden_resource_id = getattr(args, "gene_burden_resource_id", None) or _default_resource_id(args.organism, "gene_burden")
    if use_bundle and gene_burden_resource_id and ctx is not None:
        gene_burden_path = resolve_resource_path(
            ctx=ctx,
            resource_id=gene_burden_resource_id,
            resource_policy=resource_policy,
            role_label="gene_burden",
            enablement_hint=(
                "Provide --resources_dir containing the splice gene-burden table, or set --gene_burden_resource_id to a valid local resource."
            ),
        )

    alias_map = read_event_alias_tsv(alias_path) if alias_path is not None else None
    ubiquity_map = read_event_ubiquity_tsv(ubiquity_path) if ubiquity_path is not None else None
    impact_map = read_event_impact_tsv(impact_path) if impact_path is not None else None
    gene_burden_map = read_gene_burden_tsv(gene_burden_path) if gene_burden_path is not None else None

    files = [
        input_file_record(args.psi_matrix_tsv, "psi_matrix_tsv"),
        input_file_record(args.sample_metadata_tsv, "sample_metadata_tsv"),
    ]
    if args.event_metadata_tsv:
        files.append(input_file_record(args.event_metadata_tsv, "event_metadata_tsv"))
    if args.coverage_matrix_tsv:
        files.append(input_file_record(args.coverage_matrix_tsv, "coverage_matrix_tsv"))
    if alias_path is not None:
        files.append(input_file_record(alias_path, "event_alias_table"))
    if ubiquity_path is not None:
        files.append(input_file_record(ubiquity_path, "event_ubiquity_table"))
    if impact_path is not None:
        files.append(input_file_record(impact_path, "event_impact_table"))
    if gene_burden_path is not None:
        files.append(input_file_record(gene_burden_path, "gene_burden_table"))

    dataset_label = str(args.dataset_label or "").strip() or Path(args.psi_matrix_tsv).name
    signature_name = str(args.signature_name or "").strip() or Path(args.psi_matrix_tsv).stem
    cfg = SpliceMatrixWorkflowConfig(
        converter_name="splice_event_matrix",
        out_dir=out_dir,
        organism=args.organism,
        genome_build=args.genome_build,
        dataset_label=dataset_label,
        signature_name=signature_name,
        psi_matrix_tsv=args.psi_matrix_tsv,
        sample_metadata_tsv=args.sample_metadata_tsv,
        event_metadata_tsv=args.event_metadata_tsv,
        coverage_matrix_tsv=args.coverage_matrix_tsv,
        matrix_format=args.matrix_format,
        sample_id_column=args.sample_id_column,
        group_column=args.group_column,
        condition_column=args.condition_column,
        study_contrast=args.study_contrast,
        condition_a=args.condition_a,
        condition_b=args.condition_b,
        min_samples_per_condition=int(args.min_samples_per_condition),
        effect_metric=args.effect_metric,
        missing_value_policy=args.missing_value_policy,
        min_present_per_condition=int(args.min_present_per_condition),
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
        delta_psi_soft_floor=float(getattr(args, "delta_psi_soft_floor", 0.05)),
        delta_psi_soft_floor_mode=getattr(args, "delta_psi_soft_floor_mode", "auto"),
        event_dup_policy=args.event_dup_policy,
        gene_aggregation=args.gene_aggregation,
        gene_topk_events=int(args.gene_topk_events),
        gene_burden_penalty_mode=getattr(args, "gene_burden_penalty_mode", "auto"),
        min_gene_burden_penalty=float(getattr(args, "min_gene_burden_penalty", 0.35)),
        locus_density_penalty_mode=getattr(args, "locus_density_penalty_mode", "none"),
        locus_density_window_bp=int(getattr(args, "locus_density_window_bp", 20000000)),
        locus_density_top_n=int(getattr(args, "locus_density_top_n", 20)),
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
    return run_splice_event_matrix_workflow(
        cfg=cfg,
        alias_map=alias_map,
        ubiquity_map=ubiquity_map,
        impact_map=impact_map,
        gene_burden_map=gene_burden_map,
        input_files=files,
        resources_info=resources_info,
    )
