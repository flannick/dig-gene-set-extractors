from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace

from geneset_extractors.workflows.gtex_runtime_common import write_workflow_provenance_graph
from geneset_extractors.workflows.motrpac_common import prepare_tissue_inputs, read_tsv, write_json, write_tsv
from geneset_extractors.workflows.motrpac_released_dea import LEGACY_TISSUE_TERMS
from geneset_extractors.workflows.motrpac_timepoint import run as run_motrpac_timepoint
from geneset_extractors.workflows.motrpac_timewise import run as run_motrpac_timewise
from geneset_extractors.workflows.motrpac_training import run as run_motrpac_training


def _resolve_raw_counts_tsv(raw_counts_dir: Path, raw_counts_object: str) -> Path:
    object_name = str(raw_counts_object).strip()
    candidates = [
        raw_counts_dir / f"{object_name}.tsv.gz",
        raw_counts_dir / "raw_counts_by_tissue" / f"{object_name}.tsv.gz",
    ]
    for candidate in candidates:
        if candidate.exists() and candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "Missing raw counts TSV for "
        f"{object_name}. Looked in: {', '.join(str(path) for path in candidates)}"
    )


def _build_signed_term_rows_from_pooled(
    deg_rows: list[dict[str, str]],
    *,
    tissue_id: str,
    padj_max: float,
) -> list[dict[str, str]]:
    base_term = LEGACY_TISSUE_TERMS.get(tissue_id, tissue_id)
    term = f"{base_term}_Consensus"
    out_rows: list[dict[str, str]] = []
    for row in deg_rows:
        gene_id = str(row.get("gene_id", "")).strip()
        gene_symbol = str(row.get("gene_symbol", "")).strip()
        if not gene_id or not gene_symbol:
            continue
        try:
            padj = float(str(row.get("padj", "")).strip())
            logfc = float(str(row.get("logFC", "")).strip())
        except ValueError:
            continue
        if padj <= 0.0 or padj > padj_max or logfc == 0.0:
            continue
        sign = 1.0 if logfc > 0 else -1.0
        signed_score = -math.log10(padj) * sign
        out_rows.append(
            {"term": term, "gene_id": gene_id, "gene_symbol": gene_symbol, "score": str(signed_score), "sign": str(sign)}
        )
    return out_rows


def _build_signed_term_rows_from_stratified(
    deg_rows: list[dict[str, str]],
    *,
    tissue_id: str,
    padj_max: float,
) -> list[dict[str, str]]:
    base_term = LEGACY_TISSUE_TERMS.get(tissue_id, tissue_id)
    out_rows: list[dict[str, str]] = []
    for row in deg_rows:
        comparison_id = str(row.get("comparison_id", "")).strip()
        gene_id = str(row.get("gene_id", "")).strip()
        gene_symbol = str(row.get("gene_symbol", "")).strip()
        if not comparison_id or not gene_id or not gene_symbol:
            continue
        try:
            padj = float(str(row.get("padj", "")).strip())
            logfc = float(str(row.get("logFC", "")).strip())
        except ValueError:
            continue
        if padj <= 0.0 or padj > padj_max or logfc == 0.0:
            continue
        parts = comparison_id.split("_")
        if len(parts) == 3:
            sex = parts[1].capitalize()
            timepoint = parts[2].upper()
            term = f"{base_term}_{sex}_{timepoint}"
        elif len(parts) == 2:
            timepoint = parts[1].upper()
            term = f"{base_term}_{timepoint}"
        else:
            term = comparison_id
        sign = 1.0 if logfc > 0 else -1.0
        signed_score = -math.log10(padj) * sign
        out_rows.append(
            {"term": term, "gene_id": gene_id, "gene_symbol": gene_symbol, "score": str(signed_score), "sign": str(sign)}
        )
    return out_rows


def run(args) -> dict[str, object]:
    raw_counts_dir = Path(args.raw_counts_dir).resolve()
    transcript_metadata_tsv = Path(args.transcript_metadata_tsv).resolve()
    phenotype_metadata_tsv = Path(args.phenotype_metadata_tsv).resolve()
    feature_to_gene_tsv = Path(args.feature_to_gene_tsv).resolve()
    rat_to_human_tsv = Path(args.rat_to_human_tsv).resolve()
    tissue_list_tsv = Path(args.tissue_list_tsv).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    tissue_rows = read_tsv(tissue_list_tsv)
    workflow_mode = str(args.workflow_mode).strip()
    signed_term_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, object]] = []
    input_paths: list[tuple[Path, str]] = [
        (transcript_metadata_tsv, "transcript_metadata_tsv"),
        (phenotype_metadata_tsv, "phenotype_metadata_tsv"),
        (feature_to_gene_tsv, "feature_to_gene_tsv"),
        (rat_to_human_tsv, "rat_to_human_tsv"),
        (tissue_list_tsv, "tissue_list_tsv"),
    ]

    staging_root = out_dir / "staging"
    for tissue_row in tissue_rows:
        tissue_id = str(tissue_row.get("tissue_id", "")).strip()
        tissue_label = str(tissue_row.get("tissue_label", "")).strip()
        transcript_tissue_label = str(tissue_row.get("transcript_tissue_label", "")).strip()
        raw_counts_object = str(tissue_row.get("raw_counts_object", "")).strip()
        if not tissue_id or not tissue_label or not transcript_tissue_label or not raw_counts_object:
            continue
        counts_tsv = _resolve_raw_counts_tsv(raw_counts_dir, raw_counts_object)
        input_paths.append((counts_tsv, "counts_tsv"))
        prepared = prepare_tissue_inputs(
            counts_tsv=counts_tsv,
            transcript_metadata_tsv=transcript_metadata_tsv,
            phenotype_metadata_tsv=phenotype_metadata_tsv,
            feature_to_gene_tsv=feature_to_gene_tsv,
            rat_to_human_tsv=rat_to_human_tsv,
            tissue_label=tissue_label,
            transcript_tissue_label=transcript_tissue_label,
        )
        tissue_stage = staging_root / tissue_id
        tissue_stage.mkdir(parents=True, exist_ok=True)
        prepared_counts = tissue_stage / "tissue_counts.tsv"
        prepared_meta = tissue_stage / "sample_metadata.tsv"
        write_tsv(prepared_counts, prepared["counts_rows"], prepared["counts_fieldnames"])
        write_tsv(prepared_meta, prepared["sample_metadata_rows"], prepared["sample_metadata_fieldnames"])
        write_json(tissue_stage / "prepare_summary.json", prepared["summary"])

        workflow_out = tissue_stage / "workflow"
        workflow_out.mkdir(parents=True, exist_ok=True)
        if workflow_mode == "pooled":
            result = run_motrpac_training(
                SimpleNamespace(
                    counts_tsv=str(prepared_counts),
                    sample_metadata_tsv=str(prepared_meta),
                    out_dir=str(workflow_out),
                    organism="human",
                    genome_build="hg38",
                    rscript_bin=str(args.rscript_bin),
                    covariates=str(getattr(args, "covariates", "sex")),
                )
            )
            deg_rows = read_tsv(Path(result["deg_tsv_path"]))
            rows_out = _build_signed_term_rows_from_pooled(deg_rows, tissue_id=tissue_id, padj_max=float(args.padj_max))
        elif workflow_mode == "stratified_sex_timepoint":
            result = run_motrpac_timewise(
                SimpleNamespace(
                    counts_tsv=str(prepared_counts),
                    sample_metadata_tsv=str(prepared_meta),
                    out_dir=str(workflow_out),
                    organism="human",
                    genome_build="hg38",
                    min_samples_per_group=int(args.min_samples_per_group),
                )
            )
            deg_rows = read_tsv(Path(result["deg_long_path"]))
            rows_out = _build_signed_term_rows_from_stratified(deg_rows, tissue_id=tissue_id, padj_max=float(args.padj_max))
        elif workflow_mode == "stratified_timepoint":
            result = run_motrpac_timepoint(
                SimpleNamespace(
                    counts_tsv=str(prepared_counts),
                    sample_metadata_tsv=str(prepared_meta),
                    out_dir=str(workflow_out),
                    organism="human",
                    genome_build="hg38",
                    tissue_id=tissue_id,
                    rscript_bin=str(args.rscript_bin),
                    min_samples_per_group=int(args.min_samples_per_group),
                )
            )
            deg_rows = read_tsv(Path(result["deg_long_path"]))
            rows_out = _build_signed_term_rows_from_stratified(deg_rows, tissue_id=tissue_id, padj_max=float(args.padj_max))
        else:
            raise ValueError(f"Unsupported motrpac_raw_aggregated workflow_mode: {workflow_mode}")

        signed_term_rows.extend(rows_out)
        audit_rows.append(
            {
                "tissue_id": tissue_id,
                "workflow_mode": workflow_mode,
                "signed_term_rows": len(rows_out),
            }
        )

    signed_term_path = out_dir / "motrpac_signed_term_gene.tsv"
    audit_path = out_dir / "motrpac_processing_audit.tsv"
    write_tsv(signed_term_path, signed_term_rows, ["term", "gene_id", "gene_symbol", "score", "sign"])
    write_tsv(audit_path, audit_rows, ["tissue_id", "workflow_mode", "signed_term_rows"])
    write_workflow_provenance_graph(
        workflow_name="motrpac_raw_aggregated",
        module_name="geneset_extractors.workflows.motrpac_raw_aggregated",
        output_dir=out_dir,
        focus_output_path=signed_term_path,
        output_paths=[
            (signed_term_path, "signed_term_gene_tsv"),
            (audit_path, "processing_audit_tsv"),
        ],
        input_paths=input_paths,
        parameters={
            "workflow_mode": workflow_mode,
            "covariates": str(getattr(args, "covariates", "sex")),
            "padj_max": float(args.padj_max),
            "min_samples_per_group": int(args.min_samples_per_group),
            "n_signed_term_rows": len(signed_term_rows),
            "n_tissues": len(audit_rows),
        },
    )
    return {"n_rows": len(signed_term_rows), "out_dir": str(out_dir)}
