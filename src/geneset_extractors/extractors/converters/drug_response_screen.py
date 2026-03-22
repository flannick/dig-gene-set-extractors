from __future__ import annotations

import os
from pathlib import Path
import sys

from geneset_extractors.core.metadata import input_file_record
from geneset_extractors.core.provenance import activate_runtime_context
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
    filter_promiscuous_targets,
    load_compound_qc_tsv,
    load_drug_alias_map_tsv,
    load_target_ubiquity_tsv,
    load_target_aliases_tsv,
    normalize_drug_target_ids,
    summarize_targets_per_drug,
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


def _default_drug_alias_map_resource_id(organism: str) -> str | None:
    org = str(organism).strip().lower()
    if org == "human":
        return "drug_alias_map_human_v1"
    return None


def _default_blacklist_resource_id(organism: str) -> str | None:
    org = str(organism).strip().lower()
    if org == "human":
        return "drug_response_blacklist_human"
    return None


def _default_target_edges_resource_id(organism: str) -> str | None:
    org = str(organism).strip().lower()
    if org == "human":
        return "drug_target_edges_human_v1"
    return None


def _default_target_ubiquity_resource_id(organism: str) -> str | None:
    org = str(organism).strip().lower()
    if org == "human":
        return "target_ubiquity_human_v1"
    return None


def _default_compound_qc_resource_id(organism: str) -> str | None:
    org = str(organism).strip().lower()
    if org == "human":
        return "compound_qc_human_v1"
    return None


def _env_has_resources_dir() -> bool:
    return bool(os.getenv("GENESET_EXTRACTORS_RESOURCES_DIR") or os.getenv("OMICS2GENESET_RESOURCES_DIR"))


def _normalize_response_record_drug_ids(
    response_records: list[ResponseRecord],
    *,
    alias_map: dict[str, str] | None,
) -> tuple[list[ResponseRecord], dict[str, object]]:
    if not alias_map:
        return list(response_records), {"n_records": len(response_records), "n_alias_hits": 0, "n_collisions": 0}
    out: list[ResponseRecord] = []
    n_alias_hits = 0
    n_collisions = 0
    seen_pairs: set[tuple[str, str, float]] = set()
    for record in response_records:
        canonical = str(alias_map.get(record.drug_id, record.drug_id)).strip()
        if canonical != record.drug_id:
            n_alias_hits += 1
        normalized = ResponseRecord(sample_id=record.sample_id, drug_id=canonical, response=record.response)
        key = (normalized.sample_id, normalized.drug_id, float(normalized.response))
        if key in seen_pairs:
            n_collisions += 1
        seen_pairs.add(key)
        out.append(normalized)
    return out, {
        "n_records": len(response_records),
        "n_alias_hits": n_alias_hits,
        "n_collisions": n_collisions,
    }


def _resolve_program_preset(args) -> str:
    name = str(getattr(args, "program_preset", "connectable") or "connectable").strip().lower()
    if name == "default":
        return "connectable"
    if name not in {"connectable", "broad_pharmacology"}:
        raise ValueError(f"Unsupported program_preset: {name}")
    return name


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
    activate_runtime_context("drug_response_screen", getattr(args, "provenance_overlay_json", None))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    resource_policy = str(args.resource_policy)
    if resource_policy not in {"skip", "fail"}:
        raise ValueError(f"Unsupported resource_policy: {resource_policy}")

    program_preset = _resolve_program_preset(args)

    use_resources = bool(
        args.resources_dir
        or args.resources_manifest
        or args.target_aliases_resource_id
        or args.drug_blacklist_resource_id
        or getattr(args, "drug_alias_map_resource_id", None)
        or getattr(args, "target_edges_resource_id", None)
        or getattr(args, "target_ubiquity_resource_id", None)
        or getattr(args, "compound_qc_resource_id", None)
        or _env_has_resources_dir()
    )
    resource_ctx = load_resource_context(args.resources_manifest, args.resources_dir) if use_resources else None

    def _bundle_resolve(resource_id: str | None, role_label: str, enablement_hint: str, *, fallback_behavior: str) -> str | None:
        if resource_ctx is None:
            return None
        resolved = resolve_resource_path(
            ctx=resource_ctx,
            resource_id=resource_id,
            resource_policy=resource_policy,
            role_label=role_label,
            enablement_hint=enablement_hint,
        )
        if resolved is None and resource_id:
            print(
                f"warning: bundle_resource_missing resource={resource_id} action={fallback_behavior}",
                file=sys.stderr,
            )
        return str(resolved) if resolved is not None else None

    target_aliases_resource_id = args.target_aliases_resource_id or _default_alias_resource_id(args.organism)
    target_aliases_path = args.target_aliases_tsv
    if target_aliases_path is None and resource_ctx is not None:
        target_aliases_path = _bundle_resolve(
            target_aliases_resource_id,
            "target_aliases",
            "Provide --target_aliases_tsv or set --resources_dir/--resources_manifest with a valid alias resource.",
            fallback_behavior="continue_without_target_aliases",
        )

    blacklist_resource_id = args.drug_blacklist_resource_id or _default_blacklist_resource_id(args.organism)
    drug_blacklist_path = args.drug_blacklist_tsv
    if drug_blacklist_path is None and resource_ctx is not None:
        drug_blacklist_path = _bundle_resolve(
            blacklist_resource_id,
            "drug_blacklist",
            "Provide --drug_blacklist_tsv or set --resources_dir/--resources_manifest with a valid blacklist resource.",
            fallback_behavior="continue_without_blacklist",
        )

    drug_alias_map_resource_id = getattr(args, "drug_alias_map_resource_id", None) or _default_drug_alias_map_resource_id(args.organism)
    drug_alias_map_path = _bundle_resolve(
        drug_alias_map_resource_id,
        "drug_alias_map",
        "Provide a compatible --resources_dir with drug_alias_map_human_v1 or disable bundle-aware normalization.",
        fallback_behavior="continue_with_input_drug_ids",
    )
    target_edges_resource_id = getattr(args, "target_edges_resource_id", None) or _default_target_edges_resource_id(args.organism)
    target_edges_bundle_path = None
    if not args.drug_targets_tsv:
        target_edges_bundle_path = _bundle_resolve(
            target_edges_resource_id,
            "drug_target_edges",
            "Provide --drug_targets_tsv or a bundle with drug_target_edges_human_v1.",
            fallback_behavior="fall_back_to_prism_treatment_info_or_error",
        )
    target_ubiquity_resource_id = getattr(args, "target_ubiquity_resource_id", None) or _default_target_ubiquity_resource_id(args.organism)
    target_ubiquity_bundle_path = _bundle_resolve(
        target_ubiquity_resource_id,
        "target_ubiquity",
        "Provide a compatible bundle with target_ubiquity_human_v1.",
        fallback_behavior="compute_subset_derived_target_ubiquity",
    )
    use_compound_qc_bundle = bool(args.use_compound_qc_bundle)
    if getattr(args, "use_compound_qc_bundle", None) is None:
        use_compound_qc_bundle = bool(program_preset == "connectable")
    compound_qc_resource_id = getattr(args, "compound_qc_resource_id", None) or _default_compound_qc_resource_id(args.organism)
    compound_qc_bundle_path = None
    if use_compound_qc_bundle:
        compound_qc_bundle_path = _bundle_resolve(
            compound_qc_resource_id,
            "compound_qc",
            "Provide a compatible bundle with compound_qc_human_v1 or disable --use_compound_qc_bundle.",
            fallback_behavior="continue_without_compound_qc_bundle",
        )

    aliases: dict[str, str] | None = None
    aliases_summary: dict[str, object] | None = None
    if target_aliases_path:
        aliases, aliases_summary = load_target_aliases_tsv(
            target_aliases_path,
            alias_column=args.alias_column,
            gene_symbol_column=args.alias_gene_symbol_column,
            delimiter=args.alias_delimiter,
        )

    drug_alias_map: dict[str, str] | None = None
    drug_alias_summary: dict[str, object] | None = None
    if drug_alias_map_path:
        drug_alias_map, drug_alias_summary = load_drug_alias_map_tsv(drug_alias_map_path)

    target_ubiquity_bundle: dict[str, float] | None = None
    target_ubiquity_bundle_summary: dict[str, object] | None = None
    if target_ubiquity_bundle_path:
        target_ubiquity_bundle, target_ubiquity_bundle_summary = load_target_ubiquity_tsv(target_ubiquity_bundle_path)

    compound_qc_rows: dict[str, dict[str, object]] | None = None
    compound_qc_summary: dict[str, object] | None = None
    if compound_qc_bundle_path:
        compound_qc_rows, compound_qc_summary = load_compound_qc_tsv(compound_qc_bundle_path)

    used_by_id = {str(record["id"]): record for record in (resource_ctx.used if resource_ctx is not None else [])}

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
        if not args.drug_targets_tsv and not target_edges_bundle_path:
            raise ValueError("Generic mode requires --drug_targets_tsv")
        response_records, response_summary = read_response_tsv(
            path=args.response_tsv,
            sample_id_column=args.sample_id_column,
            drug_id_column=args.drug_id_column,
            response_column=args.response_column,
            delimiter=args.response_delimiter,
        )
        target_rows, target_io_summary = read_drug_targets_tsv(
            path=(args.drug_targets_tsv or target_edges_bundle_path),
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
        files.append(
            input_file_record(
                args.drug_targets_tsv or target_edges_bundle_path,
                "drug_targets_tsv",
                resource_record=None if args.drug_targets_tsv else used_by_id.get(str(target_edges_resource_id)),
            )
        )
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
        if args.drug_targets_tsv or target_edges_bundle_path:
            target_rows, target_io_summary = read_drug_targets_tsv(
                path=(args.drug_targets_tsv or target_edges_bundle_path),
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
            target_summary = {
                "build": target_build_summary,
                "fallback": "prism_treatment_info_targets",
            }
            print(
                "warning: falling back to PRISM treatment-info target text because no explicit --drug_targets_tsv "
                "or bundle drug_target_edges resource was available.",
                file=sys.stderr,
            )
        files.append(input_file_record(args.prism_matrix_csv, "prism_matrix_csv"))
        files.append(input_file_record(args.prism_treatment_info_csv, "prism_treatment_info_csv"))
        files.append(input_file_record(args.prism_cell_line_info_csv, "prism_cell_line_info_csv"))
        if args.drug_targets_tsv or target_edges_bundle_path:
            files.append(
                input_file_record(
                    args.drug_targets_tsv or target_edges_bundle_path,
                    "drug_targets_tsv",
                    resource_record=None if args.drug_targets_tsv else used_by_id.get(str(target_edges_resource_id)),
                )
            )

    alias_apply_summary: dict[str, object] | None = None
    target_id_apply_summary: dict[str, object] | None = None
    if drug_alias_map:
        response_records, alias_apply_summary = _normalize_response_record_drug_ids(
            response_records,
            alias_map=drug_alias_map,
        )
        drug_targets, target_id_apply_summary = normalize_drug_target_ids(
            drug_targets,
            alias_map=drug_alias_map,
        )

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
        files.append(
            input_file_record(
                drug_blacklist_path,
                "drug_blacklist_tsv",
                resource_record=None if args.drug_blacklist_tsv else used_by_id.get(str(blacklist_resource_id)),
            )
        )

    bundle_drug_weights: dict[str, float] = {}
    compound_qc_effect_summary: dict[str, object] | None = None
    if compound_qc_rows and use_compound_qc_bundle:
        kept_records: list[ResponseRecord] = []
        n_drop = 0
        n_warn = 0
        n_polypharm_flag = 0
        warned_drugs: list[str] = []
        drop_ids: set[str] = set()
        for drug_id, row in compound_qc_rows.items():
            recommended = str(row.get("recommended_use", "keep")).strip().lower()
            pan_toxic = bool(row.get("pan_toxic_flag", False))
            if bool(getattr(args, "include_blacklisted_compounds", False)):
                pass
            elif recommended == "drop" or pan_toxic:
                drop_ids.add(drug_id)
            elif recommended == "warn":
                n_warn += 1
                warned_drugs.append(drug_id)
            if bool(row.get("polypharm_flag", False)):
                n_polypharm_flag += 1
                bundle_drug_weights[drug_id] = float(bundle_drug_weights.get(drug_id, 1.0)) * 0.5
        if drop_ids:
            for record in response_records:
                if record.drug_id in drop_ids:
                    n_drop += 1
                    continue
                kept_records.append(record)
            response_records = kept_records
            drug_targets = {drug_id: targets for drug_id, targets in drug_targets.items() if drug_id not in drop_ids}
        for drug_id in warned_drugs[:20]:
            print(
                f"warning: compound_qc_bundle recommends caution for drug={drug_id}; recommended_use=warn",
                file=sys.stderr,
            )
        compound_qc_effect_summary = {
            "source": "bundle",
            "n_drugs_flagged_warn": n_warn,
            "n_response_rows_dropped": n_drop,
            "n_drugs_dropped": len(drop_ids),
            "n_drugs_polypharm_flag": n_polypharm_flag,
            "include_blacklisted_compounds": bool(getattr(args, "include_blacklisted_compounds", False)),
        }

    if aliases_summary is not None:
        target_summary["aliases"] = aliases_summary
        files.append(
            input_file_record(
                target_aliases_path,
                "target_aliases_tsv",
                resource_record=None if args.target_aliases_tsv else used_by_id.get(str(target_aliases_resource_id)),
            )
        )
    if drug_alias_summary is not None:
        target_summary["drug_alias_map"] = drug_alias_summary
        if alias_apply_summary is not None:
            target_summary["drug_alias_map"]["applied_to_response"] = alias_apply_summary
        if target_id_apply_summary is not None:
            target_summary["drug_alias_map"]["applied_to_targets"] = target_id_apply_summary
        files.append(
            input_file_record(
                drug_alias_map_path,
                "drug_alias_map_tsv",
                resource_record=used_by_id.get(str(drug_alias_map_resource_id)),
            )
        )
    if target_ubiquity_bundle_summary is not None:
        target_summary["target_ubiquity_bundle"] = target_ubiquity_bundle_summary
        files.append(
            input_file_record(
                target_ubiquity_bundle_path,
                "target_ubiquity_bundle_tsv",
                resource_record=used_by_id.get(str(target_ubiquity_resource_id)),
            )
        )
    if compound_qc_summary is not None:
        target_summary["compound_qc_bundle"] = compound_qc_summary
        if compound_qc_effect_summary is not None:
            target_summary["compound_qc_bundle"]["effects"] = compound_qc_effect_summary
        files.append(
            input_file_record(
                compound_qc_bundle_path,
                "compound_qc_bundle_tsv",
                resource_record=used_by_id.get(str(compound_qc_resource_id)),
            )
        )

    drug_targets, promisc_summary, promisc_warnings = filter_promiscuous_targets(
        drug_targets=drug_targets,
        max_targets_per_drug=int(args.max_targets_per_drug),
        policy=("warn" if program_preset == "broad_pharmacology" else str(args.target_promiscuity_policy)),
    )
    for msg in promisc_warnings:
        print(msg, file=sys.stderr)
    target_summary["promiscuity"] = promisc_summary
    target_summary["targets_per_drug"] = summarize_targets_per_drug(drug_targets)
    target_summary["n_drugs_zero_targets_after_cleaning"] = sum(1 for v in drug_targets.values() if not v)
    target_summary["program_preset"] = program_preset
    if target_edges_bundle_path and not args.drug_targets_tsv:
        target_summary["target_source"] = "bundle:drug_target_edges_human_v1"
    elif args.drug_targets_tsv:
        target_summary["target_source"] = "explicit:drug_targets_tsv"
    else:
        target_summary["target_source"] = "prism_treatment_info"
    if program_preset == "connectable" and target_ubiquity_bundle is None and str(args.target_ubiquity_penalty) == "idf":
        print(
            "warning: target_ubiquity_penalty=idf is using a subset-derived table because target_ubiquity_human_v1 was not available.",
            file=sys.stderr,
        )

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
    response_ubiquity_penalty = str(args.response_ubiquity_penalty or args.ubiquity_penalty)
    gmt_split_signed = args.gmt_split_signed
    if gmt_split_signed is None:
        gmt_split_signed = bool(contrast_method in {"group_vs_rest", "case_control"})

    resources_info = build_resources_info(resource_ctx) if resource_ctx is not None else None
    if resources_info is not None:
        target_summary["resources"] = resources_info
    if resource_ctx is not None:
        for record in resource_ctx.used:
            files.append(input_file_record(str(record["path"]), f"resource:{record['id']}", resource_record=record))

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
        min_group_size=int(args.min_group_size),
        scoring_model=args.scoring_model,
        sparse_alpha=float(args.sparse_alpha),
        response_ubiquity_penalty=response_ubiquity_penalty,
        target_ubiquity_penalty=str(args.target_ubiquity_penalty),
        ubiquity_tau=float(args.ubiquity_tau),
        ubiquity_epsilon=float(args.ubiquity_epsilon),
        polypharm_downweight=bool(args.polypharm_downweight),
        polypharm_t0=int(args.polypharm_t0),
        external_drug_weights=(bundle_drug_weights or None),
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
        gmt_format=str(args.gmt_format),
        emit_small_gene_sets=bool(args.emit_small_gene_sets),
        gtf=args.gtf,
        gtf_source=args.gtf_source,
        gtf_gene_id_field=args.gtf_gene_id_field,
    )

    result = run_drug_response_workflow(
        cfg=cfg,
        response_records=response_records,
        drug_targets=drug_targets,
        target_ubiquity_bundle=target_ubiquity_bundle,
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
