from __future__ import annotations

from pathlib import Path

from geneset_extractors.core.metadata import input_file_record
from geneset_extractors.extractors.rnaseq.deg_scoring import sanitize_name_component
from geneset_extractors.extractors.rnaseq.sc_program_workflow import (
    SCRNAProgramsWorkflowConfig,
    read_program_loadings,
    run_sc_programs_workflow,
)


def _resolve_source(args) -> tuple[str, str]:
    provided: list[tuple[str, str]] = []
    if args.program_loadings_tsv:
        provided.append(("program_loadings_tsv", str(args.program_loadings_tsv)))
    if args.cnmf_gene_spectra_tsv:
        provided.append(("cnmf_gene_spectra_tsv", str(args.cnmf_gene_spectra_tsv)))
    if args.schpf_gene_scores_tsv:
        provided.append(("schpf_gene_scores_tsv", str(args.schpf_gene_scores_tsv)))
    if len(provided) != 1:
        raise ValueError(
            "Provide exactly one input source: --program_loadings_tsv, --cnmf_gene_spectra_tsv, or --schpf_gene_scores_tsv"
        )
    return provided[0]


def _resolve_requested_format(args, source_kind: str, source_path: str) -> str:
    requested = str(args.loadings_format).strip()
    if source_kind == "cnmf_gene_spectra_tsv":
        if requested == "auto":
            if args.cnmf_kind in {"tpm", "score"}:
                return f"cnmf_gene_spectra_{args.cnmf_kind}"
            lower = Path(source_path).name.lower()
            if "gene_spectra_tpm" in lower:
                return "cnmf_gene_spectra_tpm"
            if "gene_spectra_score" in lower:
                return "cnmf_gene_spectra_score"
        return requested
    if source_kind == "schpf_gene_scores_tsv":
        if requested == "auto":
            return "schpf_gene_scores"
        return requested
    return requested


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    source_kind, source_path = _resolve_source(args)
    requested_format = _resolve_requested_format(args, source_kind, source_path)
    programs, parse_summary = read_program_loadings(
        path=source_path,
        loadings_format=requested_format,
        gene_id_column=args.gene_id_column,
        program_id_column=args.program_id_column,
        loading_column=args.loading_column,
        transpose=bool(args.transpose),
    )

    signature_name = str(args.signature_name or "").strip()
    if not signature_name or signature_name == "programs":
        signature_name = sanitize_name_component(Path(source_path).stem)

    dataset_label = str(args.dataset_label or "").strip()
    if not dataset_label:
        dataset_label = Path(source_path).name

    files = [input_file_record(source_path, source_kind)]
    if args.gtf:
        files.append(input_file_record(args.gtf, "gtf"))

    cfg = SCRNAProgramsWorkflowConfig(
        converter_name="rna_sc_programs",
        out_dir=out_dir,
        organism=args.organism,
        genome_build=args.genome_build,
        dataset_label=dataset_label,
        signature_name=signature_name,
        source_path=source_path,
        source_kind=source_kind,
        loadings_format=str(parse_summary.get("effective_loadings_format", requested_format)),
        score_transform=args.score_transform,
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
        exclude_gene_regex=args.exclude_gene_regex,
        disable_default_excludes=bool(args.disable_default_excludes),
        gtf=args.gtf,
        gtf_gene_id_field=args.gtf_gene_id_field,
        gtf_source=args.gtf_source,
        program_id_prefix=args.program_id_prefix,
    )
    result = run_sc_programs_workflow(
        cfg=cfg,
        programs=programs,
        parse_summary=parse_summary,
        input_files=files,
    )
    return {
        "n_groups": int(result.get("n_groups", 0)),
        "out_dir": str(out_dir),
    }
