from __future__ import annotations

from pathlib import Path

from omics2geneset.core.metadata import input_file_record
from omics2geneset.io.gtf import read_genes_from_gtf
from omics2geneset.methylation.resource_utils import (
    build_resources_info,
    load_resource_context,
    resolve_resource_path,
)
from omics2geneset.methylation.workflow import (
    MethylationWorkflowConfig,
    parse_dmr_region_features,
    read_tsv_rows,
    run_methylation_workflow,
)


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fieldnames, rows = read_tsv_rows(args.dmr_tsv)
    genes = read_genes_from_gtf(args.gtf)
    gene_chroms = {str(g.chrom) for g in genes}

    resource_policy = str(args.resource_policy)
    if resource_policy not in {"skip", "fail"}:
        raise ValueError(f"Unsupported resource_policy: {resource_policy}")

    need_resources = bool(args.enhancer_resource_id)
    ctx = load_resource_context(args.resources_manifest, args.resources_dir) if need_resources else None

    enhancer_bed = args.enhancer_bed
    if not enhancer_bed and args.enhancer_resource_id and ctx is not None:
        resolved = resolve_resource_path(
            ctx=ctx,
            resource_id=args.enhancer_resource_id,
            resource_policy=resource_policy,
            role_label="enhancer_bed",
            enablement_hint=(
                "Provide --enhancer_bed or add an enhancer BED resource and pass --enhancer_resource_id."
            ),
        )
        enhancer_bed = str(resolved) if resolved is not None else None

    features, parse_summary = parse_dmr_region_features(
        fieldnames=fieldnames,
        rows=rows,
        chrom_column=args.chrom_column,
        start_column=args.start_column,
        end_column=args.end_column,
        delta_column=args.delta_column,
        padj_column=args.padj_column,
        pvalue_column=args.pvalue_column,
        score_mode=args.score_mode,
        neglog10p_eps=float(args.neglog10p_eps),
        neglog10p_cap=float(args.neglog10p_cap),
        drop_sex_chrom=bool(args.drop_sex_chrom),
        gene_chroms=gene_chroms,
    )

    files = [
        input_file_record(args.dmr_tsv, "dmr_tsv"),
        input_file_record(args.gtf, "gtf"),
    ]
    if args.region_gene_links_tsv:
        files.append(input_file_record(args.region_gene_links_tsv, "region_gene_links_tsv"))
    if enhancer_bed:
        files.append(input_file_record(enhancer_bed, "enhancer_bed"))
    if ctx is not None:
        for record in ctx.used:
            files.append(input_file_record(str(record["path"]), f"resource:{record['id']}"))

    cfg = MethylationWorkflowConfig(
        converter_name="methylation_dmr_regions",
        out_dir=out_dir,
        organism=args.organism,
        genome_build=args.genome_build,
        dataset_label=(str(args.dataset_label).strip() if args.dataset_label else Path(args.dmr_tsv).name),
        gtf=args.gtf,
        gtf_source=args.gtf_source,
        program_preset=args.program_preset,
        program_methods=args.program_methods,
        link_method=args.link_method,
        region_gene_links_tsv=args.region_gene_links_tsv,
        promoter_upstream_bp=int(args.promoter_upstream_bp),
        promoter_downstream_bp=int(args.promoter_downstream_bp),
        max_distance_bp=(int(args.max_distance_bp) if args.max_distance_bp is not None else None),
        decay_length_bp=int(args.decay_length_bp),
        max_genes_per_peak=int(args.max_genes_per_peak),
        score_transform=args.score_transform,
        normalize=args.normalize,
        aggregation=args.aggregation,
        select=args.select,
        top_k=int(args.top_k),
        quantile=float(args.quantile),
        min_score=float(args.min_score),
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
        resource_policy=resource_policy,
        score_mode=args.score_mode,
        distal_mode=args.distal_mode,
        enhancer_bed=enhancer_bed,
    )

    resources_info = build_resources_info(ctx) if ctx is not None else None
    result = run_methylation_workflow(
        cfg=cfg,
        features=features,
        parse_summary=parse_summary,
        input_files=files,
        resources_info=resources_info,
    )
    return result
