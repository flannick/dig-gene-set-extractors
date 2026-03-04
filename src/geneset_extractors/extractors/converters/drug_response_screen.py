from __future__ import annotations

from pathlib import Path
import sys

from geneset_extractors.core.metadata import input_file_record
from geneset_extractors.extractors.drug_response.io import (
    ResponseRecord,
    read_case_control_tsv,
    read_drug_id_list_tsv,
    read_drug_targets_tsv,
    read_groups_tsv,
    read_prism_cell_line_info_csv,
    read_prism_matrix_csv,
    read_prism_treatment_info_csv,
    read_response_tsv,
    read_sample_metadata_tsv,
)
from geneset_extractors.extractors.drug_response.normalize import resolve_response_direction
from geneset_extractors.extractors.drug_response.resource_utils import (
    build_resources_info,
    load_resource_context,
    resolve_resource_path,
)
from geneset_extractors.extractors.drug_response.targets import (
    apply_drug_blacklist,
    build_drug_targets_from_table_rows,
    build_drug_targets_from_target_text,
    load_target_aliases_tsv,
)
from geneset_extractors.extractors.drug_response.workflow import (
    DrugResponseWorkflowConfig,
    run_drug_response_workflow,
)


def _load_groups_from_metadata(
    *,
    sample_metadata: dict[str, dict[str, str]] | None,
    group_column: str,
) -> dict[str, str]:
    if not sample_metadata:
        return {}
    out: dict[str, str] = {}
    for sample_id, row in sample_metadata.items():
        value = str(row.get(group_column, "")).strip()
        if value:
            out[sample_id] = value
    return out


def _resolve_contrast_method(args, groups_by_sample: dict[str, str]) -> str:
    explicit = str(args.contrast_method or "").strip()
    if explicit:
        return explicit
    if groups_by_sample:
        return "group_vs_rest"
    return "none"


def _default_alias_resource_id(organism: str) -> str | None:
    org = str(organism).strip().lower()
    if org == "human":
        return "drug_response_target_aliases_human"
    return None


def _default_blacklist_resource_id(organism: str) -> str | None:
    org = str(organism).strip().lower()
    if org == "human":
        return "drug_response_blacklist_human"
    return None


def _resolve_mode(args) -> str:
    has_generic = bool(args.response_tsv)
    has_prism = bool(args.prism_matrix_csv or args.prism_treatment_info_csv or args.prism_cell_line_info_csv)
    if has_generic and has_prism:
        raise ValueError(
            "Provide either generic mode inputs (--response_tsv + --drug_targets_tsv) "
            "or PRISM mode inputs (--prism_matrix_csv + --prism_treatment_info_csv + --prism_cell_line_info_csv), not both."
        )
    if has_generic:
        return "generic"
    if has_prism:
        return "prism"
    raise ValueError(
        "No input mode selected. Provide generic inputs (--response_tsv and --drug_targets_tsv) "
        "or PRISM inputs (--prism_matrix_csv, --prism_treatment_info_csv, --prism_cell_line_info_csv)."
    )


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    resource_policy = str(args.resource_policy)
    if resource_policy not in {"skip", "fail"}:
        raise ValueError(f"Unsupported resource_policy: {resource_policy}")

    use_resources = bool(
        args.resources_dir
        or args.resources_manifest
        or args.target_aliases_resource_id
        or args.drug_blacklist_resource_id
    )
    resource_ctx = load_resource_context(args.resources_manifest, args.resources_dir) if use_resources else None

    target_aliases_path = args.target_aliases_tsv
    if target_aliases_path is None and resource_ctx is not None:
        target_aliases_resource_id = args.target_aliases_resource_id or _default_alias_resource_id(args.organism)
        resolved = resolve_resource_path(
            ctx=resource_ctx,
            resource_id=target_aliases_resource_id,
            resource_policy=resource_policy,
            role_label="target_aliases",
            enablement_hint=(
                "Provide --target_aliases_tsv or set --resources_dir/--resources_manifest "
                "with a valid alias resource."
            ),
        )
        target_aliases_path = str(resolved) if resolved is not None else None

    drug_blacklist_path = args.drug_blacklist_tsv
    if drug_blacklist_path is None and resource_ctx is not None:
        blacklist_resource_id = args.drug_blacklist_resource_id or _default_blacklist_resource_id(args.organism)
        resolved = resolve_resource_path(
            ctx=resource_ctx,
            resource_id=blacklist_resource_id,
            resource_policy=resource_policy,
            role_label="drug_blacklist",
            enablement_hint=(
                "Provide --drug_blacklist_tsv or set --resources_dir/--resources_manifest "
                "with a valid blacklist resource."
            ),
        )
        drug_blacklist_path = str(resolved) if resolved is not None else None

    aliases: dict[str, str] | None = None
    aliases_summary: dict[str, object] | None = None
    if target_aliases_path:
        aliases, aliases_summary = load_target_aliases_tsv(
            target_aliases_path,
            alias_column=args.alias_column,
            gene_symbol_column=args.alias_gene_symbol_column,
            delimiter=args.alias_delimiter,
        )

    mode = _resolve_mode(args)
    response_records: list[ResponseRecord] = []
    response_summary: dict[str, object] = {}
    sample_metadata: dict[str, dict[str, str]] | None = None
    sample_metadata_summary: dict[str, object] | None = None
    drug_targets: dict[str, dict[str, float]] = {}
    target_summary: dict[str, object] = {}
    files = []

    if mode == "generic":
        if not args.response_tsv:
            raise ValueError("Generic mode requires --response_tsv")
        if not args.drug_targets_tsv:
            raise ValueError("Generic mode requires --drug_targets_tsv")
        response_records, response_summary = read_response_tsv(
            path=args.response_tsv,
            sample_id_column=args.sample_id_column,
            drug_id_column=args.drug_id_column,
            response_column=args.response_column,
            delimiter=args.response_delimiter,
        )
        target_rows, target_io_summary = read_drug_targets_tsv(
            path=args.drug_targets_tsv,
            drug_id_column=args.targets_drug_id_column,
            gene_symbol_column=args.targets_gene_symbol_column,
            weight_column=args.targets_weight_column,
            source_column=args.targets_source_column,
            delimiter=args.targets_delimiter,
        )
        drug_targets, target_build_summary = build_drug_targets_from_table_rows(
            target_rows,
            aliases=aliases,
        )
        target_summary = {"io": target_io_summary, "build": target_build_summary}
        files.append(input_file_record(args.response_tsv, "response_tsv"))
        files.append(input_file_record(args.drug_targets_tsv, "drug_targets_tsv"))
        if args.sample_metadata_tsv:
            sample_metadata, sample_metadata_summary = read_sample_metadata_tsv(
                path=args.sample_metadata_tsv,
                sample_id_column=args.sample_metadata_sample_id_column,
                delimiter=args.sample_metadata_delimiter,
            )
            files.append(input_file_record(args.sample_metadata_tsv, "sample_metadata_tsv"))
    else:
        if not args.prism_matrix_csv or not args.prism_treatment_info_csv or not args.prism_cell_line_info_csv:
            raise ValueError(
                "PRISM mode requires --prism_matrix_csv, --prism_treatment_info_csv, and --prism_cell_line_info_csv"
            )
        column_to_drug, drug_to_target_text, treatment_summary = read_prism_treatment_info_csv(
            path=args.prism_treatment_info_csv,
            column_name_column=args.prism_column_name_column,
            broad_id_column=args.prism_broad_id_column,
            target_column=args.prism_target_column,
        )
        response_records, prism_matrix_summary = read_prism_matrix_csv(
            path=args.prism_matrix_csv,
            column_to_drug=column_to_drug,
            sample_id_column=args.prism_matrix_sample_id_column,
        )
        sample_metadata, sample_metadata_summary = read_prism_cell_line_info_csv(
            path=args.prism_cell_line_info_csv,
            sample_id_column=args.prism_cell_line_sample_id_column,
        )
        response_summary = {
            "mode": "prism",
            "treatment_info": treatment_summary,
            "matrix": prism_matrix_summary,
            "cell_line_info": sample_metadata_summary,
        }
        if args.drug_targets_tsv:
            target_rows, target_io_summary = read_drug_targets_tsv(
                path=args.drug_targets_tsv,
                drug_id_column=args.targets_drug_id_column,
                gene_symbol_column=args.targets_gene_symbol_column,
                weight_column=args.targets_weight_column,
                source_column=args.targets_source_column,
                delimiter=args.targets_delimiter,
            )
            drug_targets, target_build_summary = build_drug_targets_from_table_rows(
                target_rows,
                aliases=aliases,
            )
            target_summary = {"io": target_io_summary, "build": target_build_summary}
        else:
            drug_targets, target_build_summary = build_drug_targets_from_target_text(
                drug_to_target_text=drug_to_target_text,
                aliases=aliases,
            )
            target_summary = {"build": target_build_summary}
        files.append(input_file_record(args.prism_matrix_csv, "prism_matrix_csv"))
        files.append(input_file_record(args.prism_treatment_info_csv, "prism_treatment_info_csv"))
        files.append(input_file_record(args.prism_cell_line_info_csv, "prism_cell_line_info_csv"))
        if args.drug_targets_tsv:
            files.append(input_file_record(args.drug_targets_tsv, "drug_targets_tsv"))

    groups_by_sample: dict[str, str] = {}
    groups_summary: dict[str, object] | None = None
    if args.groups_tsv:
        groups_by_sample, groups_summary = read_groups_tsv(
            path=args.groups_tsv,
            sample_id_column=args.groups_sample_id_column,
            group_column=args.groups_group_column,
            delimiter=args.groups_delimiter,
        )
        files.append(input_file_record(args.groups_tsv, "groups_tsv"))
    else:
        if sample_metadata:
            default_group_column = str(args.group_column or "").strip()
            if default_group_column:
                groups_by_sample = _load_groups_from_metadata(
                    sample_metadata=sample_metadata,
                    group_column=default_group_column,
                )
                groups_summary = {
                    "source": "sample_metadata",
                    "group_column": default_group_column,
                    "n_assignments": len(groups_by_sample),
                    "n_groups": len(set(groups_by_sample.values())) if groups_by_sample else 0,
                }
            if mode == "prism" and not groups_by_sample:
                groups_by_sample = _load_groups_from_metadata(
                    sample_metadata=sample_metadata,
                    group_column="primary_tissue",
                )
                groups_summary = {
                    "source": "prism_cell_line_info",
                    "group_column": "primary_tissue",
                    "n_assignments": len(groups_by_sample),
                    "n_groups": len(set(groups_by_sample.values())) if groups_by_sample else 0,
                }

    case_control_by_sample: dict[str, bool] | None = None
    case_control_summary: dict[str, object] | None = None
    if args.case_control_tsv:
        case_control_by_sample, case_control_summary = read_case_control_tsv(
            path=args.case_control_tsv,
            sample_id_column=args.case_control_sample_id_column,
            is_case_column=args.case_control_is_case_column,
            delimiter=args.case_control_delimiter,
        )
        files.append(input_file_record(args.case_control_tsv, "case_control_tsv"))

    if drug_blacklist_path:
        blacklist, blacklist_summary = read_drug_id_list_tsv(drug_blacklist_path, delimiter=args.drug_blacklist_delimiter)
        if blacklist:
            response_records = [r for r in response_records if str(r.drug_id) not in blacklist]
            drug_targets, blacklist_apply_summary = apply_drug_blacklist(
                drug_targets=drug_targets,
                blacklist=blacklist,
            )
            target_summary["blacklist"] = {
                **blacklist_summary,
                **blacklist_apply_summary,
            }
        files.append(input_file_record(drug_blacklist_path, "drug_blacklist_tsv"))

    if aliases_summary is not None:
        target_summary["aliases"] = aliases_summary
        files.append(input_file_record(target_aliases_path, "target_aliases_tsv"))

    contrast_method = _resolve_contrast_method(args, groups_by_sample)
    if contrast_method in {"group_mean", "group_vs_rest"} and not groups_by_sample:
        raise ValueError(
            f"contrast_method={contrast_method} requires groups. Provide --groups_tsv or sample metadata with --group_column."
        )
    if contrast_method == "case_control" and case_control_by_sample is None:
        raise ValueError("contrast_method=case_control requires --case_control_tsv")

    response_direction = resolve_response_direction(
        response_metric=args.response_metric,
        response_direction=args.response_direction,
    )
    gmt_split_signed = args.gmt_split_signed
    if gmt_split_signed is None:
        gmt_split_signed = bool(contrast_method in {"group_vs_rest", "case_control"})

    resources_info = build_resources_info(resource_ctx) if resource_ctx is not None else None
    if resources_info is not None:
        target_summary["resources"] = resources_info
    if resource_ctx is not None:
        for record in resource_ctx.used:
            files.append(input_file_record(str(record["path"]), f"resource:{record['id']}"))

    dataset_label = str(args.dataset_label or "").strip()
    if not dataset_label:
        dataset_label = Path(args.response_tsv or args.prism_matrix_csv).name

    cfg = DrugResponseWorkflowConfig(
        converter_name="drug_response_screen",
        out_dir=out_dir,
        organism=args.organism,
        genome_build=args.genome_build,
        dataset_label=dataset_label,
        response_metric=args.response_metric,
        response_direction=response_direction,
        response_transform=args.response_transform,
        contrast_method=contrast_method,
        case_control_within_group=bool(args.case_control_within_group),
        scoring_model=args.scoring_model,
        sparse_alpha=float(args.sparse_alpha),
        ubiquity_penalty=args.ubiquity_penalty,
        ubiquity_tau=float(args.ubiquity_tau),
        ubiquity_epsilon=float(args.ubiquity_epsilon),
        polypharm_downweight=bool(args.polypharm_downweight),
        polypharm_t0=int(args.polypharm_t0),
        max_programs=int(args.max_programs),
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
        gmt_split_signed=gmt_split_signed,
        emit_small_gene_sets=bool(args.emit_small_gene_sets),
        gtf=args.gtf,
        gtf_source=args.gtf_source,
        gtf_gene_id_field=args.gtf_gene_id_field,
    )

    result = run_drug_response_workflow(
        cfg=cfg,
        response_records=response_records,
        drug_targets=drug_targets,
        response_summary={
            **response_summary,
            "sample_metadata": sample_metadata_summary,
            "groups": groups_summary,
            "case_control": case_control_summary,
        },
        target_summary=target_summary,
        groups_by_sample=(groups_by_sample or None),
        case_control_by_sample=case_control_by_sample,
        input_files=files,
    )

    if resources_info:
        print(
            "warning: drug_response resources summary: "
            f"used={len(resources_info.get('used', []))} "
            f"missing={len(resources_info.get('missing', []))}",
            file=sys.stderr,
        )

    return {
        "n_groups": int(result.get("n_groups", 0)),
        "out_dir": str(out_dir),
    }
