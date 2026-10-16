from __future__ import annotations

from pathlib import Path

from geneset_extractors.extractors.proteomics.public_prepare import run_public_prepare
from geneset_extractors.workflows.gtex_runtime_common import write_workflow_provenance_graph

_OUTPUT_ROLES = [
    ("ptm_matrix_tsv", "ptm_matrix"),
    ("protein_matrix_tsv", "protein_matrix"),
    ("sample_metadata_tsv", "sample_metadata"),
    ("sample_id_map_tsv", "sample_id_map"),
    ("site_id_map_tsv", "site_id_map"),
    ("bundle_source_row_tsv", "bundle_source_row"),
]

_DESCRIPTION = (
    "Analysis step that standardizes CPTAC/CDAP phosphoproteome and proteome reports into "
    "PTM site and protein matrices with harmonized sample ids and emits ptm_matrix.tsv."
)


def _emit_prepare_provenance_graph(summary: dict[str, object], args) -> Path | None:
    out_dir = Path(str(summary["out_dir"]))
    outputs = summary.get("outputs", {})
    if not isinstance(outputs, dict) or not outputs.get("ptm_matrix_tsv"):
        return None

    input_paths = [
        (Path(str(rec["path"])), role)
        for role, rec in sorted((summary.get("source_files") or {}).items())
        if isinstance(rec, dict) and rec.get("path")
    ]
    output_paths = [
        (Path(str(outputs[key])), role)
        for key, role in _OUTPUT_ROLES
        if outputs.get(key)
    ]
    sample_annotations = getattr(args, "sample_annotations_tsv", None)
    if sample_annotations:
        sa_path = Path(str(sample_annotations))
        if all(existing != sa_path for existing, _role in input_paths):
            input_paths.append((sa_path, "sample_annotations"))
    focus = out_dir / "ptm_matrix.tsv"
    return write_workflow_provenance_graph(
        workflow_name="ptm_prepare_public",
        module_name="geneset_extractors.workflows.ptm_prepare_public",
        output_dir=out_dir,
        focus_output_path=focus,
        output_paths=output_paths,
        input_paths=input_paths,
        parameters={
            "ptm_type": args.ptm_type,
            "organism": args.organism,
            "study_id": summary.get("study_id"),
            "parser_profile": summary.get("parser_profile"),
        },
        description=_DESCRIPTION,
    )


def run(args) -> dict[str, object]:
    summary = run_public_prepare(
        input_mode=args.input_mode,
        ptm_report_tsv=args.ptm_report_tsv,
        protein_report_tsv=args.protein_report_tsv,
        sample_design_tsv=args.sample_design_tsv,
        sample_annotations_tsv=args.sample_annotations_tsv,
        pdc_manifest_tsv=args.pdc_manifest_tsv,
        source_dir=args.source_dir,
        out_dir=args.out_dir,
        organism=args.organism,
        ptm_type=args.ptm_type,
        study_id=args.study_id,
        study_label=args.study_label,
        assay_type_policy=args.assay_type_policy,
        min_phospho_like_fraction=args.min_phospho_like_fraction,
        max_k_fraction=args.max_k_fraction,
    )
    _emit_prepare_provenance_graph(summary, args)
    return summary
