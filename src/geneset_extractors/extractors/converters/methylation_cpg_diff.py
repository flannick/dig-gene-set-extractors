from __future__ import annotations

import os
from pathlib import Path
import sys

from geneset_extractors.core.metadata import input_file_record
from geneset_extractors.core.provenance import activate_runtime_context
from geneset_extractors.io.gtf import read_genes_from_gtf
from geneset_extractors.extractors.methylation.resource_utils import (
    build_resources_info,
    load_resource_context,
    resolve_resource_path,
)
from geneset_extractors.extractors.methylation.workflow import (
    MethylationWorkflowConfig,
    parse_cpg_diff_features,
    read_probe_blacklist,
    read_probe_manifest_tsv,
    read_tsv_rows,
    run_methylation_workflow,
)


def _default_probe_manifest_resource_id(genome_build: str, array_type: str) -> str | None:
    gb = str(genome_build).strip().lower()
    arr = str(array_type).strip().lower()
    if gb in {"hg19", "grch37"} and arr == "450k":
        return "methylation_probe_manifest_450k_hg19"
    if gb in {"hg19", "grch37"} and arr == "epic":
        return "methylation_probe_manifest_epic_hg19"
    if gb in {"hg38", "grch38"} and arr == "450k":
        return "methylation_probe_manifest_450k_hg38"
    if gb in {"hg38", "grch38"} and arr == "epic":
        return "methylation_probe_manifest_epic_hg38"
    return None


def run(args) -> dict[str, object]:
    activate_runtime_context("methylation_cpg_diff", getattr(args, "provenance_overlay_json", None))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fieldnames, rows = read_tsv_rows(args.cpg_tsv)
    genes = read_genes_from_gtf(args.gtf)
    gene_chroms = {str(g.chrom) for g in genes}

    resource_policy = str(args.resource_policy)
    if resource_policy not in {"skip", "fail"}:
        raise ValueError(f"Unsupported resource_policy: {resource_policy}")

    has_resource_location = bool(
        args.resources_dir
        or args.resources_manifest
        or os.getenv("GENESET_EXTRACTORS_RESOURCES_DIR")
        or os.getenv("OMICS2GENESET_RESOURCES_DIR")
    )
    need_resources = bool(
        has_resource_location
        or args.probe_manifest_resource_id
        or args.enhancer_resource_id
    )
    ctx = load_resource_context(args.resources_manifest, args.resources_dir) if need_resources else None

    probe_manifest_path: str | None = args.probe_manifest_tsv
    probe_manifest_source = "explicit_tsv" if probe_manifest_path else "none"
    attempted_probe_manifest_resource_id: str | None = None
    if not probe_manifest_path:
        probe_manifest_resource_id = args.probe_manifest_resource_id or _default_probe_manifest_resource_id(
            args.genome_build,
            args.array_type,
        )
        attempted_probe_manifest_resource_id = probe_manifest_resource_id
        if probe_manifest_resource_id and ctx is not None:
            resolved = resolve_resource_path(
                ctx=ctx,
                resource_id=probe_manifest_resource_id,
                resource_policy=resource_policy,
                role_label="probe_manifest",
                enablement_hint=(
                    "Provide --probe_manifest_tsv (probe_id->coordinate mapping), or point "
                    "--resources_dir to a bundle containing this resource id."
                ),
            )
            if resolved is not None:
                probe_manifest_path = str(resolved)
                probe_manifest_source = f"resource:{probe_manifest_resource_id}"

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

    probe_manifest = None
    if probe_manifest_path:
        probe_manifest = read_probe_manifest_tsv(
            probe_manifest_path,
            probe_id_column=args.manifest_probe_id_column,
            chrom_column=args.manifest_chrom_column,
            pos_column=args.manifest_pos_column,
            start_column=args.manifest_start_column,
            end_column=args.manifest_end_column,
            pos_is_0based=bool(args.manifest_pos_is_0based),
        )

    probe_blacklist = read_probe_blacklist(args.probe_blacklist_tsv)

    features, parse_summary = parse_cpg_diff_features(
        fieldnames=fieldnames,
        rows=rows,
        probe_id_column=args.probe_id_column,
        chrom_column=args.chrom_column,
        pos_column=args.pos_column,
        start_column=args.start_column,
        end_column=args.end_column,
        delta_column=args.delta_column,
        padj_column=args.padj_column,
        pvalue_column=args.pvalue_column,
        score_mode=args.score_mode,
        delta_orientation=args.delta_orientation,
        neglog10p_eps=float(args.neglog10p_eps),
        neglog10p_cap=float(args.neglog10p_cap),
        input_pos_is_0based=bool(args.input_pos_is_0based),
        probe_manifest=probe_manifest,
        probe_blacklist=probe_blacklist,
        drop_sex_chrom=bool(args.drop_sex_chrom),
        gene_chroms=gene_chroms,
    )

    unresolved = int(parse_summary.get("n_probe_unresolved", 0) or 0)
    parse_summary["probe_manifest_source"] = probe_manifest_source
    parse_summary["probe_manifest_resource_id_attempted"] = attempted_probe_manifest_resource_id
    parse_summary["probe_manifest_path"] = probe_manifest_path
    if unresolved > 0:
        default_hint = ""
        if attempted_probe_manifest_resource_id:
            default_hint = (
                f" Default auto-resolve resource id for this run: "
                f"{attempted_probe_manifest_resource_id}."
            )
        msg = (
            f"{unresolved} probe_id rows could not be resolved to coordinates. "
            "Provide --probe_manifest_tsv, or set --resources_dir so the converter can auto-resolve "
            "a probe manifest resource."
            + default_hint
        )
        if resource_policy == "fail":
            raise ValueError(msg)
        print(f"warning: {msg}", file=sys.stderr)

    files = [
        input_file_record(args.cpg_tsv, "cpg_tsv"),
        input_file_record(args.gtf, "gtf"),
    ]
    if probe_manifest_path:
        files.append(input_file_record(probe_manifest_path, "probe_manifest_tsv"))
    if args.probe_blacklist_tsv:
        files.append(input_file_record(args.probe_blacklist_tsv, "probe_blacklist_tsv"))
    if args.region_gene_links_tsv:
        files.append(input_file_record(args.region_gene_links_tsv, "region_gene_links_tsv"))
    if enhancer_bed:
        files.append(input_file_record(enhancer_bed, "enhancer_bed"))
    if ctx is not None:
        for record in ctx.used:
            files.append(input_file_record(str(record["path"]), f"resource:{record['id']}", resource_record=record))

    cfg = MethylationWorkflowConfig(
        converter_name="methylation_cpg_diff",
        out_dir=out_dir,
        organism=args.organism,
        genome_build=args.genome_build,
        dataset_label=(str(args.dataset_label).strip() if args.dataset_label else Path(args.cpg_tsv).name),
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
        delta_orientation=args.delta_orientation,
        distal_mode=args.distal_mode,
        enhancer_bed=enhancer_bed,
        exclude_gene_symbol_regex=args.exclude_gene_symbol_regex,
        exclude_gene_symbols_tsv=args.exclude_gene_symbols_tsv,
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
