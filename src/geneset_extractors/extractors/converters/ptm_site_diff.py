from __future__ import annotations

import os
from pathlib import Path

from geneset_extractors.core.metadata import input_file_record
from geneset_extractors.extractors.proteomics.ptm_site_diff_workflow import (
    PTMWorkflowConfig,
    read_site_alias_tsv,
    read_site_ubiquity_tsv,
    read_tsv_rows,
    run_ptm_site_diff_workflow,
)
from geneset_extractors.extractors.proteomics.resource_utils import (
    build_resources_info,
    load_resource_context,
    resolve_resource_path,
)


DEFAULT_ALIAS_RESOURCE_ID = "phosphosite_aliases_human_v1"
DEFAULT_UBIQUITY_RESOURCE_ID = "phosphosite_ubiquity_human_v1"



def _default_alias_resource_id(organism: str, ptm_type: str) -> str | None:
    if str(organism).strip().lower() == "human" and str(ptm_type).strip().lower() == "phospho":
        return DEFAULT_ALIAS_RESOURCE_ID
    return None



def _default_ubiquity_resource_id(organism: str, ptm_type: str) -> str | None:
    if str(organism).strip().lower() == "human" and str(ptm_type).strip().lower() == "phospho":
        return DEFAULT_UBIQUITY_RESOURCE_ID
    return None



def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fieldnames, rows = read_tsv_rows(args.ptm_tsv)
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
        use_bundle and (has_resource_location or args.site_alias_resource_id or args.site_ubiquity_resource_id)
    )
    ctx = load_resource_context(args.resources_manifest, args.resources_dir) if need_resources else None

    alias_path = None
    alias_resource_id = args.site_alias_resource_id or _default_alias_resource_id(args.organism, args.ptm_type)
    if use_bundle and alias_resource_id and ctx is not None:
        alias_path = resolve_resource_path(
            ctx=ctx,
            resource_id=alias_resource_id,
            resource_policy=resource_policy,
            role_label="site_alias",
            enablement_hint=(
                "Provide --resources_dir containing the PTM alias table, or set --site_alias_resource_id "
                "to a valid local resource."
            ),
        )
    ubiquity_path = None
    ubiquity_resource_id = args.site_ubiquity_resource_id or _default_ubiquity_resource_id(args.organism, args.ptm_type)
    if use_bundle and ubiquity_resource_id and ctx is not None:
        ubiquity_path = resolve_resource_path(
            ctx=ctx,
            resource_id=ubiquity_resource_id,
            resource_policy=resource_policy,
            role_label="site_ubiquity",
            enablement_hint=(
                "Provide --resources_dir containing the PTM ubiquity table, or set --site_ubiquity_resource_id "
                "to a valid local resource."
            ),
        )

    alias_map = read_site_alias_tsv(alias_path) if alias_path is not None else None
    ubiquity_map = read_site_ubiquity_tsv(ubiquity_path) if ubiquity_path is not None else None

    files = [input_file_record(args.ptm_tsv, "ptm_tsv")]
    if alias_path is not None:
        files.append(input_file_record(alias_path, "site_alias_table"))
    if ubiquity_path is not None:
        files.append(input_file_record(ubiquity_path, "site_ubiquity_table"))

    dataset_label = str(args.dataset_label or "").strip() or Path(args.ptm_tsv).name
    signature_name = str(args.signature_name or "").strip() or Path(args.ptm_tsv).stem
    cfg = PTMWorkflowConfig(
        converter_name="ptm_site_diff",
        out_dir=out_dir,
        organism=args.organism,
        genome_build=args.genome_build,
        dataset_label=dataset_label,
        signature_name=signature_name,
        ptm_type=args.ptm_type,
        site_id_column=args.site_id_column,
        site_group_column=args.site_group_column,
        gene_id_column=args.gene_id_column,
        gene_symbol_column=args.gene_symbol_column,
        protein_accession_column=args.protein_accession_column,
        residue_column=args.residue_column,
        position_column=args.position_column,
        score_column=args.score_column,
        stat_column=args.stat_column,
        logfc_column=args.logfc_column,
        padj_column=args.padj_column,
        pvalue_column=args.pvalue_column,
        localization_prob_column=args.localization_prob_column,
        peptide_count_column=args.peptide_count_column,
        protein_logfc_column=args.protein_logfc_column,
        protein_stat_column=args.protein_stat_column,
        score_mode=args.score_mode,
        score_transform=args.score_transform,
        protein_adjustment=args.protein_adjustment,
        protein_adjustment_lambda=float(args.protein_adjustment_lambda),
        confidence_weight_mode=args.confidence_weight_mode,
        min_localization_prob=float(args.min_localization_prob),
        site_dup_policy=args.site_dup_policy,
        gene_aggregation=args.gene_aggregation,
        gene_topk_sites=int(args.gene_topk_sites),
        ambiguous_gene_policy=args.ambiguous_gene_policy,
        select=args.select,
        top_k=int(args.top_k),
        quantile=float(args.quantile),
        min_score=float(args.min_score),
        normalize=args.normalize,
        emit_full=bool(args.emit_full),
        emit_gmt=bool(args.emit_gmt),
        gmt_out=args.gmt_out,
        gmt_prefer_symbol=bool(args.gmt_prefer_symbol),
        gmt_require_symbol=bool(args.gmt_require_symbol),
        gmt_biotype_allowlist=args.gmt_biotype_allowlist,
        gmt_min_genes=int(args.gmt_min_genes),
        gmt_max_genes=int(args.gmt_max_genes),
        gmt_topk_list=args.gmt_topk_list,
        gmt_mass_list=args.gmt_mass_list,
        gmt_split_signed=bool(args.gmt_split_signed),
        emit_small_gene_sets=bool(args.emit_small_gene_sets),
        neglog10p_cap=float(args.neglog10p_cap),
        neglog10p_eps=float(args.neglog10p_eps),
    )
    resources_info = build_resources_info(ctx) if ctx is not None else None
    return run_ptm_site_diff_workflow(
        cfg=cfg,
        fieldnames=fieldnames,
        rows=rows,
        alias_map=alias_map,
        ubiquity_map=ubiquity_map,
        input_files=files,
        resources_info=resources_info,
    )
