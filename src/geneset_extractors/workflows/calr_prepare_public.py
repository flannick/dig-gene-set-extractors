from __future__ import annotations

from geneset_extractors.extractors.calorimetry.public_prepare import run_public_prepare


def run(args) -> dict[str, object]:
    return run_public_prepare(
        studies_tsv=args.studies_tsv,
        out_dir=args.out_dir,
        organism=args.organism,
        output_gene_species=getattr(args, "output_gene_species", "human"),
        ortholog_policy=getattr(args, "ortholog_policy", "unique_only"),
        mouse_human_orthologs_tsv=getattr(args, "mouse_human_orthologs_tsv", None),
        bundle_id=args.bundle_id,
        build_bundle=args.build_bundle,
        write_distribution_artifact=args.write_distribution_artifact,
        term_templates_tsv=args.term_templates_tsv,
        phenotype_gene_edges_tsv=args.phenotype_gene_edges_tsv,
        term_hierarchy_tsv=args.term_hierarchy_tsv,
        include_packaged_term_hierarchy=args.include_packaged_term_hierarchy,
        exploratory_without_session=args.exploratory_without_session,
        min_group_size=args.min_group_size,
        mass_covariate=args.mass_covariate,
        analysis_start_hour=args.analysis_start_hour,
        analysis_end_hour=args.analysis_end_hour,
        photoperiod_lights_on_hour=args.photoperiod_lights_on_hour,
        photoperiod_hours_light=args.photoperiod_hours_light,
    )
