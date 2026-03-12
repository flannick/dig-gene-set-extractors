from __future__ import annotations

from geneset_extractors.extractors.proteomics.public_prepare import run_public_prepare


def run(args) -> dict[str, object]:
    return run_public_prepare(
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
    )

