from __future__ import annotations

from pathlib import Path
import sys

from geneset_extractors.core.metadata import input_file_record
from geneset_extractors.extractors.morphology.bundle import bundle_resources_info, load_bundle_context, resolve_bundle_manifest
from geneset_extractors.extractors.morphology.io import aggregate_profiles, read_compound_targets_auto, read_feature_schema_tsv, read_feature_stats_tsv, read_metadata_tsv, read_profiles_table
from geneset_extractors.extractors.morphology.workflow import MorphologyWorkflowConfig, run_morphology_workflow


REQUIRED_BUNDLE_KEYS = {
    "reference_profiles": "reference_profiles_tsv",
    "reference_metadata": "reference_metadata_tsv",
    "compound_targets": "compound_targets_tsv",
    "feature_schema": "feature_schema_tsv",
    "feature_stats": "feature_stats_tsv",
}


def _resolve_reference_inputs(args) -> tuple[dict[str, str | None], dict[str, object] | None, dict[str, object] | None]:
    explicit = {
        "reference_profiles": args.reference_profiles_tsv or args.reference_profiles_parquet,
        "reference_metadata": args.reference_metadata_tsv,
        "compound_targets": args.compound_targets_tsv,
        "feature_schema": args.feature_schema_tsv,
        "feature_stats": args.feature_stats_tsv,
    }
    if explicit["reference_profiles"] and explicit["reference_metadata"] and explicit["compound_targets"]:
        return explicit, None, None
    if not args.reference_bundle_id:
        raise ValueError(
            "Provide explicit reference files (--reference_profiles_tsv/--reference_profiles_parquet, --reference_metadata_tsv, --compound_targets_tsv) "
            "or set --reference_bundle_id with --resources_dir."
        )
    ctx = load_bundle_context(args.resources_manifest, args.resources_dir)
    bundle_manifest_path, bundle_manifest = resolve_bundle_manifest(
        ctx=ctx,
        bundle_id=args.reference_bundle_id,
        resource_policy=str(args.resource_policy),
    )
    if bundle_manifest is None or bundle_manifest_path is None:
        raise ValueError("Could not resolve morphology reference bundle.")
    bundle_root = bundle_manifest_path.parent
    files = bundle_manifest.get("files", {}) if isinstance(bundle_manifest, dict) else {}
    if not isinstance(files, dict):
        raise ValueError("bundle.json must contain a files object")
    resolved = dict(explicit)
    for key in REQUIRED_BUNDLE_KEYS:
        if resolved.get(key):
            continue
        rel = str(files.get(key, "")).strip()
        if not rel:
            if key in {"feature_schema", "feature_stats"}:
                continue
            raise ValueError(f"bundle.json missing required file entry: {key}")
        resolved[key] = str(bundle_root / rel)
    return resolved, bundle_manifest, bundle_resources_info(ctx)


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    reference_inputs, bundle_manifest, resources_info = _resolve_reference_inputs(args)

    query_profiles, query_summary = read_profiles_table(
        args.query_profiles_tsv,
        id_column=args.query_id_column,
        delimiter=args.query_profiles_delimiter,
    )
    query_metadata = None
    query_metadata_summary = None
    if args.query_metadata_tsv:
        query_metadata, query_metadata_summary = read_metadata_tsv(
            args.query_metadata_tsv,
            id_column=args.query_metadata_id_column,
            delimiter=args.query_metadata_delimiter,
        )
    aggregated_queries, aggregation_summary, query_membership = aggregate_profiles(
        query_profiles,
        query_metadata,
        group_by=args.group_query_by,
        aggregate=args.query_aggregate,
    )
    reference_profiles, reference_profiles_summary = read_profiles_table(
        reference_inputs["reference_profiles"],
        id_column=args.reference_id_column,
        delimiter=args.reference_profiles_delimiter,
    )
    reference_metadata, reference_metadata_summary = read_metadata_tsv(
        reference_inputs["reference_metadata"],
        id_column=args.reference_metadata_id_column,
        delimiter=args.reference_metadata_delimiter,
    )
    compound_targets, compound_targets_summary = read_compound_targets_auto(
        reference_inputs["compound_targets"],
        compound_id_column=args.compound_id_column,
        gene_symbol_column=args.compound_target_gene_symbol_column,
        weight_column=args.compound_target_weight_column,
        delimiter=args.compound_targets_delimiter,
    )
    feature_schema = None
    feature_schema_summary = None
    if reference_inputs.get("feature_schema"):
        feature_schema, feature_schema_summary = read_feature_schema_tsv(reference_inputs["feature_schema"])
    feature_stats = None
    feature_stats_summary = None
    if reference_inputs.get("feature_stats"):
        feature_stats, feature_stats_summary = read_feature_stats_tsv(reference_inputs["feature_stats"])

    files = [input_file_record(args.query_profiles_tsv, "query_profiles_tsv")]
    if args.query_metadata_tsv:
        files.append(input_file_record(args.query_metadata_tsv, "query_metadata_tsv"))
    for key, role in (("reference_profiles", "reference_profiles"), ("reference_metadata", "reference_metadata"), ("compound_targets", "compound_targets"), ("feature_schema", "feature_schema"), ("feature_stats", "feature_stats")):
        if reference_inputs.get(key):
            files.append(input_file_record(reference_inputs[key], role))

    dataset_label = str(args.dataset_label or "").strip() or Path(args.query_profiles_tsv).name
    signature_name = str(args.signature_name or "").strip() or Path(args.query_profiles_tsv).stem
    cfg = MorphologyWorkflowConfig(
        converter_name="morphology_profile_query",
        out_dir=out_dir,
        organism=args.organism,
        genome_build=args.genome_build,
        dataset_label=dataset_label,
        signature_name=signature_name,
        similarity_metric=args.similarity_metric,
        similarity_power=float(args.similarity_power),
        polarity=args.polarity,
        query_modality_column=str(args.query_modality_column),
        same_modality_first=bool(args.same_modality_first),
        cross_modality_penalty=float(args.cross_modality_penalty),
        max_reference_neighbors=int(args.max_reference_neighbors),
        adaptive_neighbors=bool(args.adaptive_neighbors),
        min_effective_neighbors=int(args.min_effective_neighbors),
        neighbor_evidence_drop_ratio=float(args.neighbor_evidence_drop_ratio),
        mutual_neighbor_filter=bool(args.mutual_neighbor_filter),
        min_similarity=float(args.min_similarity),
        control_calibration=str(args.control_calibration),
        hubness_penalty=str(args.hubness_penalty),
        gene_recurrence_penalty=str(args.gene_recurrence_penalty),
        min_specificity_confidence_to_emit_opposite=str(args.min_specificity_confidence_to_emit_opposite),
        compound_weight=float(args.compound_weight),
        genetic_weight=float(args.genetic_weight),
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
        gmt_format=str(args.gmt_format),
        emit_small_gene_sets=bool(args.emit_small_gene_sets),
    )
    parse_summary = {
        "query_profiles": query_summary,
        "query_metadata": query_metadata_summary,
        "query_aggregation": aggregation_summary,
        "reference_profiles": reference_profiles_summary,
        "reference_metadata": reference_metadata_summary,
        "compound_targets": compound_targets_summary,
        "feature_schema": feature_schema_summary,
        "feature_stats": feature_stats_summary,
        "bundle_manifest": bundle_manifest,
    }
    if query_metadata is not None and isinstance(bundle_manifest, dict):
        bundle_contexts = bundle_manifest.get("contexts", {})
        if isinstance(bundle_contexts, dict):
            bundle_modalities = {str(key).strip().lower() for key in bundle_contexts}
            query_modalities = {
                str(row.get("perturbation_type", "")).strip().lower()
                for row in query_metadata.values()
                if str(row.get("perturbation_type", "")).strip()
            }
            missing_modalities = sorted(query_modalities - bundle_modalities)
            if missing_modalities:
                print(
                    f"warning: query metadata contains modalities absent from bundle: {missing_modalities}",
                    file=sys.stderr,
                )
    result = run_morphology_workflow(
        cfg=cfg,
        query_profiles=aggregated_queries,
        reference_profiles=reference_profiles,
        reference_metadata=reference_metadata,
        compound_targets=compound_targets,
        feature_schema=feature_schema,
        feature_stats=feature_stats,
        query_membership=query_membership,
        parse_summary=parse_summary,
        input_files=files,
        resources_info=resources_info,
        exclude_query_ids_from_reference=bool(args.exclude_query_ids_from_reference),
        query_metadata_rows=query_metadata,
    )
    return {"n_groups": int(result.get("n_groups", 0)), "out_dir": str(out_dir)}
