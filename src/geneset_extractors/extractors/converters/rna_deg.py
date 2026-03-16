from __future__ import annotations

from pathlib import Path

from geneset_extractors.core.metadata import input_file_record
from geneset_extractors.core.provenance import activate_runtime_context
from geneset_extractors.extractors.rnaseq.deg_scoring import read_deg_tsv, sanitize_name_component
from geneset_extractors.extractors.rnaseq.deg_workflow import DEGWorkflowConfig, run_deg_workflow


def run(args) -> dict[str, object]:
    activate_runtime_context("rna_deg", getattr(args, "provenance_overlay_json", None))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fieldnames, rows = read_deg_tsv(args.deg_tsv)
    signature_name = str(args.signature_name or "").strip()
    if not signature_name or signature_name == "contrast":
        signature_name = sanitize_name_component(Path(args.deg_tsv).stem)

    files = [input_file_record(args.deg_tsv, "deg_tsv")]
    if args.gtf:
        files.append(input_file_record(args.gtf, "gtf"))

    cfg = DEGWorkflowConfig(
        converter_name="rna_deg",
        out_dir=out_dir,
        organism=args.organism,
        genome_build=args.genome_build,
        signature_name=signature_name,
        deg_tsv_label=Path(args.deg_tsv).name,
        comparison_label=None,
        gene_id_column=args.gene_id_column,
        gene_symbol_column=args.gene_symbol_column,
        stat_column=args.stat_column,
        logfc_column=args.logfc_column,
        padj_column=args.padj_column,
        pvalue_column=args.pvalue_column,
        score_column=args.score_column,
        score_mode=args.score_mode,
        neglog10p_cap=args.neglog10p_cap,
        neglog10p_eps=args.neglog10p_eps,
        duplicate_gene_policy=args.duplicate_gene_policy,
        exclude_gene_regex=args.exclude_gene_regex,
        disable_default_excludes=args.disable_default_excludes,
        gtf=args.gtf,
        gtf_gene_id_field=args.gtf_gene_id_field,
        gtf_source=args.gtf_source,
        select=args.select,
        top_k=args.top_k,
        quantile=args.quantile,
        min_score=args.min_score,
        normalize=args.normalize,
        emit_full=args.emit_full,
        emit_gmt=args.emit_gmt,
        gmt_out=args.gmt_out,
        gmt_prefer_symbol=args.gmt_prefer_symbol,
        gmt_require_symbol=args.gmt_require_symbol,
        gmt_biotype_allowlist=args.gmt_biotype_allowlist,
        gmt_min_genes=args.gmt_min_genes,
        gmt_max_genes=args.gmt_max_genes,
        gmt_topk_list=args.gmt_topk_list,
        gmt_mass_list=args.gmt_mass_list,
        gmt_split_signed=args.gmt_split_signed,
        gmt_emit_abs=args.gmt_emit_abs,
        gmt_source=args.gmt_source,
        emit_small_gene_sets=args.emit_small_gene_sets,
        warn_biotype_missing=True,
    )
    result = run_deg_workflow(
        cfg=cfg,
        fieldnames=fieldnames,
        rows=rows,
        input_files=files,
    )
    return {
        "n_input_features": result["n_input_features"],
        "n_genes": result["n_genes_selected"],
        "out_dir": str(out_dir),
        "resolved_score_mode": result["resolved_score_mode"],
    }
