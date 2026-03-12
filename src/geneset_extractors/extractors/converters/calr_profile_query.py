from __future__ import annotations

from pathlib import Path

from geneset_extractors.extractors.calorimetry.bundle import bundle_file_path, bundle_resources_info, load_bundle_context, resolve_bundle_manifest
from geneset_extractors.extractors.calorimetry.workflow import CalRWorkflowConfig, run_calr_profile_workflow


REQUIRED_KEYS = {
    "reference_profiles": "reference_profiles_tsv",
    "reference_metadata": "reference_metadata_tsv",
}


def _resolve_inputs(args):
    explicit = {
        "reference_profiles_tsv": str(args.reference_profiles_tsv or "").strip() or None,
        "reference_metadata_tsv": str(args.reference_metadata_tsv or "").strip() or None,
        "feature_schema_tsv": str(args.feature_schema_tsv or "").strip() or None,
        "feature_stats_tsv": str(args.feature_stats_tsv or "").strip() or None,
    }
    if explicit["reference_profiles_tsv"] and explicit["reference_metadata_tsv"]:
        return explicit, None
    if not args.reference_bundle_id:
        raise ValueError(
            "Provide --reference_profiles_tsv and --reference_metadata_tsv, or set --reference_bundle_id with --resources_dir."
        )
    ctx = load_bundle_context(args.resources_manifest, args.resources_dir)
    bundle_manifest_path, bundle_manifest = resolve_bundle_manifest(
        ctx=ctx,
        bundle_id=str(args.reference_bundle_id),
        resource_policy=str(args.resource_policy),
    )
    if bundle_manifest is None or bundle_manifest_path is None:
        raise ValueError("Could not resolve calorimetry reference bundle.")
    resolved = dict(explicit)
    for key, output_key in REQUIRED_KEYS.items():
        if resolved[output_key]:
            continue
        resolved[output_key] = bundle_file_path(bundle_manifest_path, bundle_manifest, key)
        if not resolved[output_key]:
            raise ValueError(f"Calorimetry bundle missing required file entry: {key}")
    if not resolved["feature_schema_tsv"]:
        resolved["feature_schema_tsv"] = bundle_file_path(bundle_manifest_path, bundle_manifest, "feature_schema")
    if not resolved["feature_stats_tsv"]:
        resolved["feature_stats_tsv"] = bundle_file_path(bundle_manifest_path, bundle_manifest, "feature_stats")
    return resolved, bundle_resources_info(ctx)


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dataset_label = str(args.dataset_label or "").strip() or Path(args.calr_data_csv).name
    signature_name = str(args.signature_name or "").strip() or Path(args.calr_data_csv).stem
    cfg = CalRWorkflowConfig(
        converter_name="calr_profile_query",
        out_dir=out_dir,
        organism=args.organism,
        genome_build=args.genome_build,
        dataset_label=dataset_label,
        signature_name=signature_name,
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
        gmt_format=args.gmt_format,
        emit_small_gene_sets=bool(args.emit_small_gene_sets),
        analysis_start_hour=args.analysis_start_hour,
        analysis_end_hour=args.analysis_end_hour,
        photoperiod_lights_on_hour=args.photoperiod_lights_on_hour,
        photoperiod_hours_light=float(args.photoperiod_hours_light),
        exploratory_without_session=bool(args.exploratory_without_session),
        exclusions_tsv=args.exclusions_tsv,
        mass_covariate=args.mass_covariate,
        min_group_size=int(args.min_group_size),
        similarity_metric=args.similarity_metric,
        similarity_floor=float(args.similarity_floor),
        similarity_power=float(args.similarity_power),
        hubness_penalty=args.hubness_penalty,
        provenance_mismatch_penalty=float(args.provenance_mismatch_penalty),
        reference_bundle_id=args.reference_bundle_id,
    )
    resolved_inputs, resources_info = _resolve_inputs(args)
    return run_calr_profile_workflow(
        cfg=cfg,
        calr_data_csv=args.calr_data_csv,
        session_csv=args.session_csv,
        exclusions_tsv=args.exclusions_tsv,
        reference_profiles_tsv=resolved_inputs["reference_profiles_tsv"],
        reference_metadata_tsv=resolved_inputs["reference_metadata_tsv"],
        feature_schema_tsv=resolved_inputs["feature_schema_tsv"],
        feature_stats_tsv=resolved_inputs["feature_stats_tsv"],
        resources_info=resources_info,
    )
