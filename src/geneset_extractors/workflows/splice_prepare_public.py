from __future__ import annotations

from geneset_extractors.extractors.splicing.public_prepare import run_public_prepare


def run(args) -> dict[str, object]:
    return run_public_prepare(
        input_mode=args.input_mode,
        psi_tsv=args.psi_tsv,
        sample_annotations_tsv=args.sample_annotations_tsv,
        out_dir=args.out_dir,
        organism=args.organism,
        genome_build=args.genome_build,
        study_id=args.study_id,
        study_label=args.study_label,
        event_metadata_tsv=args.event_metadata_tsv,
        sample_id_map_tsv=args.sample_id_map_tsv,
        event_type_allowlist=args.event_type_allowlist,
        missing_psi_policy=args.missing_psi_policy,
    )
