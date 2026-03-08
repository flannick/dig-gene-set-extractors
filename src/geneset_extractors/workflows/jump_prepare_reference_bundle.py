from __future__ import annotations

import csv
import json
from pathlib import Path
from statistics import mean, median
import sys

from geneset_extractors.extractors.morphology.io import read_compound_targets_auto, read_metadata_tsv, read_profiles_table, read_target_annotations_tsv, write_bundle_manifest
from geneset_extractors.extractors.morphology.similarity import cosine_similarity


def _safe_value(value: str) -> str:
    return str(value).replace("/", "_").replace(" ", "_")


def _aggregate_consensus(rows: list[dict[str, float]], *, aggregate: str) -> dict[str, float]:
    features = sorted({feat for row in rows for feat in row})
    out: dict[str, float] = {}
    for feat in features:
        vals = [row[feat] for row in rows if feat in row]
        if not vals:
            continue
        out[feat] = float(median(vals) if aggregate == "median" else mean(vals))
    return out


def _hub_scores(consensus_profiles: dict[str, dict[str, float]], *, k: int) -> dict[str, float]:
    ids = sorted(consensus_profiles)
    features = sorted({feat for row in consensus_profiles.values() for feat in row})
    feature_centers: dict[str, float] = {}
    feature_scales: dict[str, float] = {}
    for feat in features:
        vals = [float(consensus_profiles[perturbation_id].get(feat, 0.0)) for perturbation_id in ids]
        ctr = float(mean(vals)) if vals else 0.0
        scale = float(max(1e-8, (sum((v - ctr) ** 2 for v in vals) / float(len(vals) or 1)) ** 0.5))
        feature_centers[feat] = ctr
        feature_scales[feat] = scale
    vectors = {}
    for perturbation_id in ids:
        vectors[perturbation_id] = [
            (float(consensus_profiles[perturbation_id].get(feat, 0.0)) - feature_centers[feat]) / feature_scales[feat]
            for feat in features
        ]
    out: dict[str, float] = {}
    for perturbation_id in ids:
        sims: list[float] = []
        for other_id in ids:
            if other_id == perturbation_id:
                continue
            sim = cosine_similarity(vectors[perturbation_id], vectors[other_id])
            if sim > 0.0:
                sims.append(float(sim))
        sims.sort(reverse=True)
        kept = sims[: max(1, int(k))] if sims else []
        out[perturbation_id] = float(sum(kept) / float(len(kept))) if kept else 0.0
    return out


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    profile_paths = [Path(p) for p in str(args.profile_paths).split(",") if str(p).strip()]
    if not profile_paths:
        raise ValueError("Provide at least one profile file via --profile_paths")
    experimental_metadata, exp_summary = read_metadata_tsv(args.experimental_metadata_tsv, id_column=args.metadata_join_id_column, delimiter=args.metadata_delimiter)
    compound_targets, targets_summary = read_compound_targets_auto(
        args.compound_targets_tsv,
        compound_id_column=args.compound_id_column,
        gene_symbol_column=args.compound_target_gene_symbol_column,
        weight_column=args.compound_target_weight_column,
        delimiter=args.compound_targets_delimiter,
    )
    if bool(targets_summary.get("fallback_used")):
        print("warning: compound target parsing used fallback public-JUMP target interpretation", file=sys.stderr)
    if int(targets_summary.get("n_zero_target_compounds", 0)) > 0:
        print(
            f"warning: {targets_summary['n_zero_target_compounds']} compounds had zero parsed targets in target metadata",
            file=sys.stderr,
        )
    target_annotations: dict[str, dict[str, str]] = {}
    target_annotations_summary: dict[str, object] | None = None
    if getattr(args, "target_annotations_tsv", None):
        target_annotations, target_annotations_summary = read_target_annotations_tsv(
            args.target_annotations_tsv,
            gene_symbol_column=str(getattr(args, "target_annotation_gene_symbol_column", "gene_symbol")),
            delimiter=str(getattr(args, "target_annotations_delimiter", "\t")),
        )
    all_profiles: dict[str, dict[str, float]] = {}
    n_rows = 0
    for profile_path in profile_paths:
        profiles, _summary = read_profiles_table(profile_path, id_column=args.profile_id_column, delimiter=args.profile_delimiter)
        all_profiles.update(profiles)
        n_rows += len(profiles)

    modality_timepoint_counts: dict[str, dict[str, int]] = {}
    filtered_ids: list[str] = []
    rows_by_modality_timepoint: dict[tuple[str, str], list[str]] = {}
    for profile_id in all_profiles:
        meta = experimental_metadata.get(profile_id, {})
        cell_type = str(meta.get(args.cell_type_column, "")).strip()
        timepoint = str(meta.get(args.timepoint_column, "")).strip()
        modality = str(meta.get(args.perturbation_type_column, "")).strip().lower()
        if args.cell_type_filter and cell_type != str(args.cell_type_filter):
            continue
        rows_by_modality_timepoint.setdefault((modality, timepoint), []).append(profile_id)
        modality_timepoint_counts.setdefault(modality, {})
        modality_timepoint_counts[modality][timepoint] = int(modality_timepoint_counts[modality].get(timepoint, 0)) + 1

    valid_modalities = [m for m in sorted(modality_timepoint_counts) if m in {"compound", "orf", "crispr"}]
    requested_timepoint = str(args.timepoint_filter or "").strip()
    chosen_timepoint = requested_timepoint
    if bool(args.allow_mixed_timepoints):
        chosen_contexts = {
            modality: (requested_timepoint if requested_timepoint in modality_timepoint_counts.get(modality, {}) else max(modality_timepoint_counts.get(modality, {}), key=lambda tp: modality_timepoint_counts[modality][tp], default=""))
            for modality in valid_modalities
        }
    else:
        if not chosen_timepoint:
            candidate_timepoints = sorted({tp for counts in modality_timepoint_counts.values() for tp in counts})
            chosen_timepoint = max(
                candidate_timepoints,
                key=lambda tp: (
                    sum(1 for modality in valid_modalities if tp in modality_timepoint_counts.get(modality, {})),
                    sum(int(modality_timepoint_counts.get(modality, {}).get(tp, 0)) for modality in valid_modalities),
                ),
                default="",
            )
        chosen_contexts = {modality: chosen_timepoint for modality in valid_modalities}

    included_modalities: list[str] = []
    missing_modalities: list[str] = []
    for modality in valid_modalities:
        context_timepoint = chosen_contexts.get(modality, "")
        if not context_timepoint or context_timepoint not in modality_timepoint_counts.get(modality, {}):
            missing_modalities.append(modality)
            continue
        included_modalities.append(modality)

    if bool(args.require_same_timepoint_across_modalities) and not bool(args.allow_mixed_timepoints) and chosen_timepoint:
        if missing_modalities:
            print(
                "warning: requested coherent same-timepoint morphology bundle is missing modalities at "
                f"timepoint={chosen_timepoint}; included={included_modalities} missing={missing_modalities}",
                file=sys.stderr,
            )
    if missing_modalities and not bool(args.allow_missing_modalities):
        raise ValueError(
            f"Missing required modalities for morphology bundle: included={included_modalities} missing={missing_modalities}"
        )

    for profile_id in all_profiles:
        meta = experimental_metadata.get(profile_id, {})
        if args.cell_type_filter and str(meta.get(args.cell_type_column, "")).strip() != str(args.cell_type_filter):
            continue
        modality = str(meta.get(args.perturbation_type_column, "")).strip().lower()
        timepoint = str(meta.get(args.timepoint_column, "")).strip()
        if modality not in valid_modalities:
            continue
        chosen = chosen_contexts.get(modality, "")
        if chosen and timepoint != chosen:
            continue
        filtered_ids.append(profile_id)
    grouped: dict[str, list[dict[str, float]]] = {}
    metadata_rows: dict[str, dict[str, str]] = {}
    for profile_id in filtered_ids:
        meta = experimental_metadata.get(profile_id, {})
        perturbation_id = str(meta.get(args.perturbation_id_column, profile_id)).strip()
        if not perturbation_id:
            continue
        grouped.setdefault(perturbation_id, []).append(all_profiles[profile_id])
        if perturbation_id not in metadata_rows:
            gene_symbol = str(meta.get(args.gene_symbol_column, "")).strip().upper()
            compound_id = str(meta.get(args.compound_id_column, "") or perturbation_id).strip()
            annotation_genes: set[str] = set()
            if gene_symbol:
                annotation_genes.add(gene_symbol)
            for target_gene in compound_targets.get(compound_id, {}):
                annotation_genes.add(str(target_gene).strip().upper())

            def _collapse_annotation(field: str) -> str:
                values = sorted(
                    {
                        str(target_annotations.get(gene, {}).get(field, "")).strip()
                        for gene in annotation_genes
                        if str(target_annotations.get(gene, {}).get(field, "")).strip()
                    }
                )
                return ";".join(values)

            metadata_rows[perturbation_id] = {
                "perturbation_id": perturbation_id,
                "perturbation_type": str(meta.get(args.perturbation_type_column, "")).strip(),
                "cell_type_or_line": str(meta.get(args.cell_type_column, "")).strip(),
                "timepoint": str(meta.get(args.timepoint_column, "")).strip(),
                "gene_symbol": str(meta.get(args.gene_symbol_column, "")).strip(),
                "compound_id": compound_id,
                "is_control": str(meta.get(args.is_control_column, "")).strip(),
                "control_type": str(meta.get(args.control_type_column, "")).strip(),
                "target_family": _collapse_annotation("target_family"),
                "target_class": _collapse_annotation("target_class"),
                "mechanism_label": _collapse_annotation("mechanism_label"),
                "pathway_seed": _collapse_annotation("pathway_seed"),
            }

    consensus_profiles = {perturbation_id: _aggregate_consensus(rows, aggregate=args.consensus_aggregate) for perturbation_id, rows in grouped.items()}
    feature_names = sorted({feat for prof in consensus_profiles.values() for feat in prof})
    feature_stats = []
    for feat in feature_names:
        vals = [consensus_profiles[p][feat] for p in consensus_profiles if feat in consensus_profiles[p]]
        ctr = float(mean(vals)) if vals else 0.0
        scale = float(max(1e-8, (sum((v - ctr) ** 2 for v in vals) / float(len(vals) or 1)) ** 0.5))
        feature_stats.append((feat, ctr, scale))

    for perturbation_id, rows in grouped.items():
        metadata_rows[perturbation_id]["replicate_count"] = str(len(rows))
        metadata_rows[perturbation_id]["qc_weight"] = str(min(1.0, max(0.1, float(len(rows)) / 3.0)))
    hub_scores = _hub_scores(consensus_profiles, k=int(args.hubness_k))
    for perturbation_id, hub_score in hub_scores.items():
        metadata_rows[perturbation_id]["hub_score"] = str(hub_score)

    profiles_path = out_dir / "reference_profiles.tsv.gz"
    metadata_path = out_dir / "reference_metadata.tsv.gz"
    targets_path = out_dir / "compound_targets.tsv.gz"
    target_annotations_path = out_dir / "target_annotations.tsv.gz"
    feature_schema_path = out_dir / "feature_schema.tsv.gz"
    feature_stats_path = out_dir / "feature_stats.tsv.gz"

    import gzip
    with gzip.open(profiles_path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["perturbation_id", *feature_names])
        for perturbation_id in sorted(consensus_profiles):
            writer.writerow([perturbation_id, *[consensus_profiles[perturbation_id].get(feat, 0.0) for feat in feature_names]])
    with gzip.open(metadata_path, "wt", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "perturbation_id",
            "perturbation_type",
            "cell_type_or_line",
            "timepoint",
            "gene_symbol",
            "compound_id",
            "is_control",
            "control_type",
            "qc_weight",
            "replicate_count",
            "hub_score",
            "target_family",
            "target_class",
            "mechanism_label",
            "pathway_seed",
        ]
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for perturbation_id in sorted(metadata_rows):
            writer.writerow(metadata_rows[perturbation_id])
    with gzip.open(targets_path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["compound_id", "gene_symbol", "weight"])
        for compound_id in sorted(compound_targets):
            for gene_symbol, weight in sorted(compound_targets[compound_id].items()):
                writer.writerow([compound_id, gene_symbol, weight])
    if target_annotations:
        with gzip.open(target_annotations_path, "wt", encoding="utf-8", newline="") as fh:
            annotation_fields = sorted({field for row in target_annotations.values() for field in row})
            writer = csv.DictWriter(fh, delimiter="\t", fieldnames=["gene_symbol", *annotation_fields])
            writer.writeheader()
            for gene_symbol in sorted(target_annotations):
                row = {"gene_symbol": gene_symbol}
                row.update(target_annotations[gene_symbol])
                writer.writerow(row)
    with gzip.open(feature_schema_path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["feature"])
        for feat in feature_names:
            writer.writerow([feat])
    with gzip.open(feature_stats_path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["feature", "center", "scale"])
        for feat, ctr, scale in feature_stats:
            writer.writerow([feat, ctr, scale])

    bundle_id = str(args.bundle_id or f"morphology_jump_target_pilot_{_safe_value(str(args.cell_type_filter or 'all'))}_{_safe_value(str(args.timepoint_filter or 'all'))}_v1")
    bundle_manifest = {
        "bundle_id": bundle_id,
        "version": "v1",
        "files": {
            "reference_profiles": profiles_path.name,
            "reference_metadata": metadata_path.name,
            "compound_targets": targets_path.name,
            **({"target_annotations": target_annotations_path.name} if target_annotations else {}),
            "feature_schema": feature_schema_path.name,
            "feature_stats": feature_stats_path.name,
        },
        "filters": {
            "cell_type_filter": args.cell_type_filter,
            "timepoint_filter": args.timepoint_filter,
            "profile_kind": args.profile_kind,
            "require_same_timepoint_across_modalities": bool(args.require_same_timepoint_across_modalities),
            "allow_missing_modalities": bool(args.allow_missing_modalities),
            "allow_mixed_timepoints": bool(args.allow_mixed_timepoints),
            "hubness_k": int(args.hubness_k),
        },
        "contexts": {
            modality: {
                "cell_type_or_line": args.cell_type_filter,
                "timepoint": chosen_contexts.get(modality, ""),
                "n_profiles": int(modality_timepoint_counts.get(modality, {}).get(chosen_contexts.get(modality, ""), 0)),
                "n_consensus_profiles": sum(
                    1
                    for perturbation_id, row in metadata_rows.items()
                    if str(row.get("perturbation_type", "")).strip().lower() == modality
                ),
            }
            for modality in valid_modalities
            if chosen_contexts.get(modality, "")
        },
        "summary": {
            "n_input_profiles": n_rows,
            "n_consensus_profiles": len(consensus_profiles),
            "n_features": len(feature_names),
            "experimental_metadata": exp_summary,
            "compound_targets": targets_summary,
            "target_annotations": target_annotations_summary,
            "included_modalities": included_modalities,
            "missing_modalities": missing_modalities,
            "annotation_coverage": {
                "n_annotated_genes": len(target_annotations),
                "n_refs_with_any_annotation": sum(
                    1
                    for row in metadata_rows.values()
                    if any(str(row.get(field, "")).strip() for field in ("target_family", "target_class", "mechanism_label", "pathway_seed"))
                ),
            },
            "hub_score_summary": {
                "min": min(hub_scores.values()) if hub_scores else None,
                "median": median(hub_scores.values()) if hub_scores else None,
                "max": max(hub_scores.values()) if hub_scores else None,
            },
        },
    }
    bundle_manifest_path = out_dir / f"{bundle_id}.bundle.json"
    write_bundle_manifest(bundle_manifest_path, bundle_manifest)
    bundle_summary = {
        "bundle_id": bundle_id,
        "cell_type_or_line": args.cell_type_filter,
        "timepoint": args.timepoint_filter,
        "n_consensus_profiles_total": len(consensus_profiles),
        "included_modalities": included_modalities,
        "missing_modalities": missing_modalities,
        "contexts": bundle_manifest["contexts"],
        "blocked_modalities": [
            {
                "modality": modality,
                "reason": f"no matching profiles in selected context at timepoint={chosen_contexts.get(modality, chosen_timepoint)}",
            }
            for modality in missing_modalities
        ],
    }
    (out_dir / "bundle_summary.json").write_text(json.dumps(bundle_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [
        f"bundle_id: {bundle_id}",
        f"cell_type_or_line: {args.cell_type_filter}",
        f"timepoint: {args.timepoint_filter}",
        f"n_consensus_profiles_total: {len(consensus_profiles)}",
        f"included_modalities: {included_modalities}",
        f"missing_modalities: {missing_modalities}",
        f"annotation_coverage: {bundle_manifest['summary']['annotation_coverage']}",
    ]
    for modality, context in bundle_manifest["contexts"].items():
        lines.append(
            f"modality={modality} included=yes timepoint={context['timepoint']} "
            f"n_profiles_raw={context['n_profiles']} n_consensus_profiles={context['n_consensus_profiles']}"
        )
    for row in bundle_summary["blocked_modalities"]:
        lines.append(f"blocked_modality={row['modality']} reason={row['reason']}")
    (out_dir / "bundle_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {"bundle_id": bundle_id, "bundle_manifest": str(bundle_manifest_path), "n_profiles": len(consensus_profiles)}
