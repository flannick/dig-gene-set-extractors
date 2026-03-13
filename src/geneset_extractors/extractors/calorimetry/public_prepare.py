from __future__ import annotations

import csv
import json
from pathlib import Path
from statistics import mean
from types import SimpleNamespace

from geneset_extractors.extractors.calorimetry.contrasts import build_group_contrasts
from geneset_extractors.extractors.calorimetry.features import CalorimetryFeatureConfig, extract_subject_features
from geneset_extractors.extractors.calorimetry.io import read_calr_data_csv, read_exclusions_tsv, read_session_csv
from geneset_extractors.workflows.calr_prepare_reference_bundle import run as run_prepare_bundle


def _read_tsv(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    p = Path(path)
    with p.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"TSV has no header: {p}")
        fieldnames = [str(name) for name in reader.fieldnames]
        rows = [{str(k): str(v or "") for k, v in row.items()} for row in reader]
    return fieldnames, rows


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _parse_float(value: object) -> float | None:
    text = _clean(value)
    if not text or text.lower() in {"na", "nan", "null", "none"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _resolve_path(base_dir: Path, value: str) -> str | None:
    text = _clean(value)
    if not text:
        return None
    path = Path(text)
    if not path.is_absolute():
        path = base_dir / path
    return str(path)


def _first_nonempty(row: dict[str, str], candidates: tuple[str, ...]) -> str:
    for key in candidates:
        value = _clean(row.get(key))
        if value:
            return value
    return ""


def _select_contrast(contrasts, study_row: dict[str, str]):
    requested = _first_nonempty(study_row, ("reference_group", "target_group", "perturbed_group", "contrast_id", "group"))
    if requested:
        for contrast in contrasts:
            if contrast.contrast_id == requested or contrast.group_label == requested or contrast.contrast_label == requested:
                return contrast
        raise ValueError(f"Requested reference_group/contrast {requested!r} not found among {[c.contrast_id for c in contrasts]}")
    if len(contrasts) == 1:
        return contrasts[0]
    raise ValueError(
        "Study row did not specify reference_group/contrast_id and multiple contrasts were available: "
        + ", ".join(contrast.contrast_id for contrast in contrasts)
    )


def _derive_feature_schema_and_stats(profiles: list[dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    feature_names = sorted({key for row in profiles for key in row.keys() if key != "reference_id"})
    schema_rows = [{"feature_name": feature_name} for feature_name in feature_names]
    stats_rows: list[dict[str, object]] = []
    for feature_name in feature_names:
        values = []
        for row in profiles:
            value = _parse_float(row.get(feature_name))
            if value is not None:
                values.append(float(value))
        center = float(mean(values)) if values else 0.0
        scale = float(max(1e-8, (sum((value - center) ** 2 for value in values) / float(len(values) or 1)) ** 0.5))
        stats_rows.append({"feature_name": feature_name, "center": center, "scale": scale})
    return schema_rows, stats_rows


def _qc_weight(extracted_summary: dict[str, object], contrast_metadata: dict[str, object], study_row: dict[str, str]) -> float:
    explicit = _parse_float(study_row.get("qc_weight"))
    if explicit is not None:
        return max(0.05, min(1.0, float(explicit)))
    weight = 1.0
    if str(extracted_summary.get("session_mode", "")) != "explicit":
        weight *= 0.7
    if str(extracted_summary.get("analysis_window_source", "")) == "full_trace":
        weight *= 0.8
    if int(contrast_metadata.get("n_subjects_group", 0) or 0) < 3:
        weight *= 0.85
    return round(max(0.05, min(1.0, weight)), 4)


def run_public_prepare(
    *,
    studies_tsv: str,
    out_dir: str,
    organism: str,
    bundle_id: str | None,
    build_bundle: bool,
    write_distribution_artifact: bool,
    term_templates_tsv: str | None,
    phenotype_gene_edges_tsv: str | None,
    term_hierarchy_tsv: str | None,
    include_packaged_term_hierarchy: bool,
    exploratory_without_session: bool,
    min_group_size: int,
    mass_covariate: str | None,
    analysis_start_hour: float | None,
    analysis_end_hour: float | None,
    photoperiod_lights_on_hour: float | None,
    photoperiod_hours_light: float,
) -> dict[str, object]:
    studies_path = Path(studies_tsv)
    studies_fieldnames, study_rows = _read_tsv(studies_path)
    if not study_rows:
        raise ValueError("studies_tsv contained no rows")

    out_root = Path(out_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    studies_base = studies_path.parent

    reference_profiles: list[dict[str, object]] = []
    reference_metadata: list[dict[str, object]] = []
    resolved_studies: list[dict[str, object]] = []
    all_warnings: list[str] = []

    for study_row in study_rows:
        study_id = _first_nonempty(study_row, ("study_id",))
        if not study_id:
            raise ValueError("Each studies_tsv row requires study_id")
        study_label = _first_nonempty(study_row, ("study_label",)) or study_id
        calr_data_csv = _resolve_path(studies_base, _first_nonempty(study_row, ("calr_data_csv",)))
        if not calr_data_csv:
            raise ValueError(f"Study {study_id} is missing calr_data_csv")
        session_csv = _resolve_path(studies_base, _first_nonempty(study_row, ("session_csv",)))
        exclusions_tsv = _resolve_path(studies_base, _first_nonempty(study_row, ("exclusions_tsv",)))

        gene_symbol = _first_nonempty(study_row, ("gene_symbol", "target_gene_symbol", "perturbed_gene_symbol"))
        gene_id = _first_nonempty(study_row, ("gene_id", "target_gene_id", "perturbed_gene_id")) or gene_symbol
        if not gene_symbol:
            raise ValueError(f"Study {study_id} requires gene_symbol in studies_tsv for a gene-labeled reference bundle")

        data_fieldnames, data_rows = read_calr_data_csv(calr_data_csv)
        session_fieldnames, session_rows = read_session_csv(session_csv) if session_csv else ([], None)
        exclusion_rows = read_exclusions_tsv(exclusions_tsv) if exclusions_tsv else []
        extracted = extract_subject_features(
            data_rows=data_rows,
            data_fieldnames=data_fieldnames,
            session_rows=session_rows,
            session_fieldnames=session_fieldnames,
            exclusion_rows=exclusion_rows,
            cfg=CalorimetryFeatureConfig(
                analysis_start_hour=analysis_start_hour,
                analysis_end_hour=analysis_end_hour,
                photoperiod_lights_on_hour=photoperiod_lights_on_hour,
                photoperiod_hours_light=photoperiod_hours_light,
                exclusions_tsv=exclusions_tsv,
                exploratory_without_session=exploratory_without_session,
            ),
        )
        contrasts, contrast_warnings = build_group_contrasts(
            subjects_by_id=extracted.subjects,
            explicit_mass_covariate=mass_covariate,
            min_group_size=min_group_size,
        )
        contrast = _select_contrast(contrasts, study_row)

        reference_id = _first_nonempty(study_row, ("reference_id",)) or f"{study_id}__group={contrast.contrast_id}__gene={gene_symbol}"
        profile_row: dict[str, object] = {"reference_id": reference_id}
        for feature_name, score in sorted(contrast.feature_scores.items()):
            profile_row[feature_name] = float(score)
        reference_profiles.append(profile_row)

        meta_row: dict[str, object] = {
            "reference_id": reference_id,
            "gene_id": gene_id,
            "gene_symbol": gene_symbol,
            "study_id": study_id,
            "study_label": study_label,
            "contrast_id": contrast.contrast_id,
            "contrast_label": contrast.contrast_label,
            "group_label": contrast.group_label,
            "perturbation_type": _first_nonempty(study_row, ("perturbation_type", "reference_type", "model_type")) or "genetic",
            "site": _first_nonempty(study_row, ("site", "location")),
            "sex": _first_nonempty(study_row, ("sex",)),
            "strain": _first_nonempty(study_row, ("strain",)),
            "diet": _first_nonempty(study_row, ("diet",)),
            "ambient_temperature": _first_nonempty(study_row, ("ambient_temperature", "room_temperature")),
            "acclimation_state": _first_nonempty(study_row, ("acclimation_state",)) or (
                "post_acclimation" if str(extracted.metadata_summary.get("analysis_window_source", "")) in {"session", "explicit"} else "unknown"
            ),
            "session_mode": extracted.metadata_summary.get("session_mode", ""),
            "session_group_layout_mode": extracted.metadata_summary.get("session_group_layout_mode", ""),
            "session_window_layout_mode": extracted.metadata_summary.get("session_window_layout_mode", ""),
            "analysis_window_source": extracted.metadata_summary.get("analysis_window_source", ""),
            "analysis_start": extracted.metadata_summary.get("analysis_start", ""),
            "analysis_end": extracted.metadata_summary.get("analysis_end", ""),
            "photoperiod_source": extracted.metadata_summary.get("photoperiod_source", ""),
            "mass_covariate": contrast.metadata.get("mass_covariate", ""),
            "mass_adjusted": contrast.metadata.get("mass_adjusted", False),
            "design_class": contrast.metadata.get("design_class", ""),
            "run_ids": ";".join(str(x) for x in contrast.metadata.get("run_ids", [])),
            "n_subjects_total": contrast.metadata.get("n_subjects_total", 0),
            "n_subjects_group": contrast.metadata.get("n_subjects_group", 0),
            "n_subjects_other": contrast.metadata.get("n_subjects_other", 0),
            "interaction_tested": contrast.metadata.get("interaction_tested", 0),
            "interaction_retained": contrast.metadata.get("interaction_retained", 0),
            "qc_weight": _qc_weight(extracted.metadata_summary, contrast.metadata, study_row),
            "source_calr_data_csv": calr_data_csv,
            "source_session_csv": session_csv or "",
            "source_exclusions_tsv": exclusions_tsv or "",
        }
        reference_metadata.append(meta_row)
        warnings = list(extracted.warnings) + list(contrast_warnings)
        all_warnings.extend(f"{study_id}: {warning}" for warning in warnings)
        resolved_studies.append(
            {
                "study_id": study_id,
                "study_label": study_label,
                "reference_id": reference_id,
                "gene_symbol": gene_symbol,
                "contrast_id": contrast.contrast_id,
                "session_mode": extracted.metadata_summary.get("session_mode", ""),
                "analysis_window_source": extracted.metadata_summary.get("analysis_window_source", ""),
                "session_group_layout_mode": extracted.metadata_summary.get("session_group_layout_mode", ""),
                "session_window_layout_mode": extracted.metadata_summary.get("session_window_layout_mode", ""),
                "n_subjects_total": contrast.metadata.get("n_subjects_total", 0),
                "n_subjects_group": contrast.metadata.get("n_subjects_group", 0),
                "n_features": len(contrast.feature_scores),
                "warnings_count": len(warnings),
            }
        )

    feature_schema_rows, feature_stats_rows = _derive_feature_schema_and_stats(reference_profiles)
    feature_names = [str(row["feature_name"]) for row in feature_schema_rows]

    profiles_path = out_root / "reference_profiles.tsv"
    metadata_path = out_root / "reference_metadata.tsv"
    schema_path = out_root / "feature_schema.tsv"
    stats_path = out_root / "feature_stats.tsv"
    resolved_path = out_root / "resolved_studies.tsv"
    summary_path = out_root / "prepare_summary.json"

    _write_tsv(profiles_path, ["reference_id", *feature_names], reference_profiles)
    _write_tsv(
        metadata_path,
        [
            "reference_id",
            "gene_id",
            "gene_symbol",
            "study_id",
            "study_label",
            "contrast_id",
            "contrast_label",
            "group_label",
            "perturbation_type",
            "site",
            "sex",
            "strain",
            "diet",
            "ambient_temperature",
            "acclimation_state",
            "session_mode",
            "session_group_layout_mode",
            "session_window_layout_mode",
            "analysis_window_source",
            "analysis_start",
            "analysis_end",
            "photoperiod_source",
            "mass_covariate",
            "mass_adjusted",
            "design_class",
            "run_ids",
            "n_subjects_total",
            "n_subjects_group",
            "n_subjects_other",
            "interaction_tested",
            "interaction_retained",
            "qc_weight",
            "source_calr_data_csv",
            "source_session_csv",
            "source_exclusions_tsv",
        ],
        reference_metadata,
    )
    _write_tsv(schema_path, ["feature_name"], feature_schema_rows)
    _write_tsv(stats_path, ["feature_name", "center", "scale"], feature_stats_rows)
    _write_tsv(
        resolved_path,
        [
            "study_id",
            "study_label",
            "reference_id",
            "gene_symbol",
            "contrast_id",
            "session_mode",
            "analysis_window_source",
            "session_group_layout_mode",
            "session_window_layout_mode",
            "n_subjects_total",
            "n_subjects_group",
            "n_features",
            "warnings_count",
        ],
        resolved_studies,
    )

    summary = {
        "workflow": "calr_prepare_public",
        "organism": organism,
        "studies_tsv": str(studies_path),
        "n_studies": len(study_rows),
        "n_reference_profiles": len(reference_profiles),
        "n_features": len(feature_names),
        "n_explicit_session_studies": sum(1 for row in resolved_studies if row["session_mode"] == "explicit"),
        "n_exploratory_studies": sum(1 for row in resolved_studies if row["session_mode"] != "explicit"),
        "output_files": {
            "reference_profiles_tsv": str(profiles_path),
            "reference_metadata_tsv": str(metadata_path),
            "feature_schema_tsv": str(schema_path),
            "feature_stats_tsv": str(stats_path),
            "resolved_studies_tsv": str(resolved_path),
        },
        "warnings": all_warnings,
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    result = {
        "out_dir": str(out_root),
        "n_profiles": len(reference_profiles),
        "n_features": len(feature_names),
        "prepare_summary_json": str(summary_path),
    }

    if build_bundle:
        bundle_dir = out_root / "bundle"
        bundle_args = SimpleNamespace(
            reference_profiles_tsv=str(profiles_path),
            reference_metadata_tsv=str(metadata_path),
            feature_schema_tsv=str(schema_path),
            feature_stats_tsv=str(stats_path),
            term_templates_tsv=term_templates_tsv,
            phenotype_gene_edges_tsv=phenotype_gene_edges_tsv,
            term_hierarchy_tsv=term_hierarchy_tsv,
            include_packaged_term_hierarchy=include_packaged_term_hierarchy,
            out_dir=str(bundle_dir),
            organism=organism,
            bundle_id=bundle_id or "calorimetry_public_bundle_v1",
            write_distribution_artifact=write_distribution_artifact,
            distribution_dir=None,
        )
        bundle_result = run_prepare_bundle(bundle_args)
        result.update(
            {
                "bundle_out_dir": str(bundle_dir),
                "bundle_id": bundle_result.get("bundle_id"),
                "bundle_manifest": bundle_result.get("bundle_manifest"),
            }
        )

    return result
