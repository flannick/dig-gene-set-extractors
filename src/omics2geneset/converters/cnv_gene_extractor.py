from __future__ import annotations

from pathlib import Path

from omics2geneset.cnv.seg_io import parse_segments_tsv, read_purity_tsv
from omics2geneset.cnv.seg_workflow import CNVWorkflowConfig, run_cnv_workflow
from omics2geneset.core.metadata import input_file_record


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    segments, parse_summary = parse_segments_tsv(
        path=args.segments_tsv,
        chrom_column=args.chrom_column,
        start_column=args.start_column,
        end_column=args.end_column,
        amplitude_column=args.amplitude_column,
        sample_id_column=args.sample_id_column,
        coord_system=args.coord_system,
        chrom_prefix_mode=args.chrom_prefix_mode,
    )

    purity_by_sample = None
    purity_summary = None
    if args.purity_tsv:
        purity_by_sample, purity_summary = read_purity_tsv(
            path=args.purity_tsv,
            sample_id_column=args.purity_sample_id_column,
            purity_value_column=args.purity_value_column,
        )

    files = [
        input_file_record(args.segments_tsv, "segments_tsv"),
        input_file_record(args.gtf, "gtf"),
    ]
    if args.purity_tsv:
        files.append(input_file_record(args.purity_tsv, "purity_tsv"))

    dataset_label = str(args.dataset_label or "").strip() or Path(args.segments_tsv).name
    cfg = CNVWorkflowConfig(
        converter_name="cnv_gene_extractor",
        out_dir=out_dir,
        organism=args.organism,
        genome_build=args.genome_build,
        gtf=args.gtf,
        gtf_source=args.gtf_source,
        gtf_gene_id_field=args.gtf_gene_id_field,
        dataset_label=dataset_label,
        program_preset=args.program_preset,
        program_methods=args.program_methods,
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
        coord_system=args.coord_system,
        chrom_prefix_mode=args.chrom_prefix_mode,
        min_abs_amplitude=float(args.min_abs_amplitude),
        focal_length_scale_bp=float(args.focal_length_scale_bp),
        focal_length_alpha=float(args.focal_length_alpha),
        gene_count_penalty_mode=args.gene_count_penalty,
        aggregation=args.aggregation,
        use_purity_correction=bool(args.use_purity_correction),
        purity_floor=float(args.purity_floor),
        max_abs_amplitude=float(args.max_abs_amplitude),
        emit_cohort_sets=bool(args.emit_cohort_sets),
        cohort_score_threshold=float(args.cohort_score_threshold),
        cohort_min_fraction=float(args.cohort_min_fraction),
        cohort_min_samples=int(args.cohort_min_samples),
    )

    result = run_cnv_workflow(
        cfg=cfg,
        segments=segments,
        parse_summary=parse_summary,
        purity_by_sample=purity_by_sample,
        purity_summary=purity_summary,
        input_files=files,
    )
    return result

