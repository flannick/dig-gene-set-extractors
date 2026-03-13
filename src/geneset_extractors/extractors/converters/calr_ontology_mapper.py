from __future__ import annotations

from contextlib import ExitStack
from importlib.resources import as_file
from pathlib import Path

from geneset_extractors.extractors.calorimetry.bundle import bundle_file_path, bundle_resources_info, load_bundle_context, resolve_bundle_manifest
from geneset_extractors.extractors.calorimetry.ontology import (
    packaged_resource_path,
)
from geneset_extractors.extractors.calorimetry.workflow import CalRWorkflowConfig, run_calr_ontology_workflow


DEFAULT_TERM_TEMPLATES = "calr_term_templates_mouse_v1.tsv"
DEFAULT_GENE_EDGES = "calr_phenotype_gene_edges_mouse_v1.tsv"
DEFAULT_TERM_HIERARCHY = "calr_term_hierarchy_mouse_v1.tsv"


def _resolve_inputs(args, stack: ExitStack):
    explicit_templates = str(args.term_templates_tsv or "").strip() or None
    explicit_edges = str(args.phenotype_gene_edges_tsv or "").strip() or None
    explicit_hierarchy = str(args.term_hierarchy_tsv or "").strip() or None
    explicit_orthologs = str(getattr(args, "mouse_human_orthologs_tsv", None) or "").strip() or None
    if explicit_templates and explicit_edges:
        return {
            "term_templates_tsv": explicit_templates,
            "phenotype_gene_edges_tsv": explicit_edges,
            "term_hierarchy_tsv": explicit_hierarchy,
            "mouse_human_orthologs_tsv": explicit_orthologs,
        }, None

    if args.reference_bundle_id:
        ctx = load_bundle_context(args.resources_manifest, args.resources_dir)
        bundle_manifest_path, bundle_manifest = resolve_bundle_manifest(
            ctx=ctx,
            bundle_id=str(args.reference_bundle_id),
            resource_policy=str(args.resource_policy),
        )
        if bundle_manifest is None or bundle_manifest_path is None:
            raise ValueError("Could not resolve calorimetry reference bundle.")
        resolved = {
            "term_templates_tsv": explicit_templates or bundle_file_path(bundle_manifest_path, bundle_manifest, "term_templates"),
            "phenotype_gene_edges_tsv": explicit_edges or bundle_file_path(bundle_manifest_path, bundle_manifest, "phenotype_gene_edges"),
            "term_hierarchy_tsv": explicit_hierarchy or bundle_file_path(bundle_manifest_path, bundle_manifest, "term_hierarchy"),
            "mouse_human_orthologs_tsv": explicit_orthologs or bundle_file_path(bundle_manifest_path, bundle_manifest, "mouse_human_orthologs"),
        }
        if not resolved["term_templates_tsv"] or not resolved["phenotype_gene_edges_tsv"]:
            raise ValueError("Calorimetry bundle is missing term_templates or phenotype_gene_edges entries.")
        return resolved, bundle_resources_info(ctx)

    term_templates_path = stack.enter_context(as_file(packaged_resource_path(DEFAULT_TERM_TEMPLATES)))
    phenotype_gene_edges_path = stack.enter_context(as_file(packaged_resource_path(DEFAULT_GENE_EDGES)))
    term_hierarchy_path = stack.enter_context(as_file(packaged_resource_path(DEFAULT_TERM_HIERARCHY)))
    return {
        "term_templates_tsv": str(term_templates_path),
        "phenotype_gene_edges_tsv": str(phenotype_gene_edges_path),
        "term_hierarchy_tsv": str(term_hierarchy_path),
        "mouse_human_orthologs_tsv": explicit_orthologs,
    }, {
        "manifest": "packaged_defaults",
        "resources_dir": None,
        "used": [
            {"id": "calorimetry_term_templates_mouse_v1", "method": "packaged_default", "path": str(term_templates_path)},
            {"id": "calorimetry_phenotype_gene_edges_mouse_v1", "method": "packaged_default", "path": str(phenotype_gene_edges_path)},
            {"id": "calorimetry_term_hierarchy_mouse_v1", "method": "packaged_default", "path": str(term_hierarchy_path)},
        ],
        "missing": [],
        "warnings": [],
    }


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dataset_label = str(args.dataset_label or "").strip() or Path(args.calr_data_csv).name
    signature_name = str(args.signature_name or "").strip() or Path(args.calr_data_csv).stem
    cfg = CalRWorkflowConfig(
        converter_name="calr_ontology_mapper",
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
        output_gene_species=getattr(args, "output_gene_species", "human"),
        ortholog_policy=getattr(args, "ortholog_policy", "unique_only"),
        reference_bundle_id=args.reference_bundle_id,
    )
    with ExitStack() as stack:
        resolved_inputs, resources_info = _resolve_inputs(args, stack)
        return run_calr_ontology_workflow(
            cfg=cfg,
            calr_data_csv=args.calr_data_csv,
            session_csv=args.session_csv,
            exclusions_tsv=args.exclusions_tsv,
            term_templates_tsv=resolved_inputs["term_templates_tsv"],
            phenotype_gene_edges_tsv=resolved_inputs["phenotype_gene_edges_tsv"],
            term_hierarchy_tsv=resolved_inputs["term_hierarchy_tsv"],
            mouse_human_orthologs_tsv=resolved_inputs["mouse_human_orthologs_tsv"],
            resources_info=resources_info,
        )
