from __future__ import annotations

import csv
from pathlib import Path
from statistics import mean, median

from geneset_extractors.extractors.morphology.io import read_compound_targets_tsv, read_metadata_tsv, read_profiles_table, write_bundle_manifest


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


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    profile_paths = [Path(p) for p in str(args.profile_paths).split(",") if str(p).strip()]
    if not profile_paths:
        raise ValueError("Provide at least one profile file via --profile_paths")
    experimental_metadata, exp_summary = read_metadata_tsv(args.experimental_metadata_tsv, id_column=args.metadata_join_id_column, delimiter=args.metadata_delimiter)
    compound_targets, targets_summary = read_compound_targets_tsv(
        args.compound_targets_tsv,
        compound_id_column=args.compound_id_column,
        gene_symbol_column=args.compound_target_gene_symbol_column,
        weight_column=args.compound_target_weight_column,
        delimiter=args.compound_targets_delimiter,
    )
    all_profiles: dict[str, dict[str, float]] = {}
    n_rows = 0
    for profile_path in profile_paths:
        profiles, _summary = read_profiles_table(profile_path, id_column=args.profile_id_column, delimiter=args.profile_delimiter)
        all_profiles.update(profiles)
        n_rows += len(profiles)

    filtered_ids = []
    for profile_id in all_profiles:
        meta = experimental_metadata.get(profile_id, {})
        if args.cell_type_filter and str(meta.get(args.cell_type_column, "")).strip() != str(args.cell_type_filter):
            continue
        if args.timepoint_filter and str(meta.get(args.timepoint_column, "")).strip() != str(args.timepoint_filter):
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
            metadata_rows[perturbation_id] = {
                "perturbation_id": perturbation_id,
                "perturbation_type": str(meta.get(args.perturbation_type_column, "")).strip(),
                "cell_type_or_line": str(meta.get(args.cell_type_column, "")).strip(),
                "timepoint": str(meta.get(args.timepoint_column, "")).strip(),
                "gene_symbol": str(meta.get(args.gene_symbol_column, "")).strip(),
                "compound_id": str(meta.get(args.compound_id_column, "") or perturbation_id).strip(),
                "is_control": str(meta.get(args.is_control_column, "")).strip(),
                "control_type": str(meta.get(args.control_type_column, "")).strip(),
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

    profiles_path = out_dir / "reference_profiles.tsv.gz"
    metadata_path = out_dir / "reference_metadata.tsv.gz"
    targets_path = out_dir / "compound_targets.tsv.gz"
    feature_schema_path = out_dir / "feature_schema.tsv.gz"
    feature_stats_path = out_dir / "feature_stats.tsv.gz"

    import gzip
    with gzip.open(profiles_path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["perturbation_id", *feature_names])
        for perturbation_id in sorted(consensus_profiles):
            writer.writerow([perturbation_id, *[consensus_profiles[perturbation_id].get(feat, 0.0) for feat in feature_names]])
    with gzip.open(metadata_path, "wt", encoding="utf-8", newline="") as fh:
        fieldnames = ["perturbation_id", "perturbation_type", "cell_type_or_line", "timepoint", "gene_symbol", "compound_id", "is_control", "control_type", "qc_weight", "replicate_count"]
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
            "feature_schema": feature_schema_path.name,
            "feature_stats": feature_stats_path.name,
        },
        "filters": {"cell_type_filter": args.cell_type_filter, "timepoint_filter": args.timepoint_filter, "profile_kind": args.profile_kind},
        "summary": {"n_input_profiles": n_rows, "n_consensus_profiles": len(consensus_profiles), "n_features": len(feature_names), "experimental_metadata": exp_summary, "compound_targets": targets_summary},
    }
    bundle_manifest_path = out_dir / f"{bundle_id}.bundle.json"
    write_bundle_manifest(bundle_manifest_path, bundle_manifest)
    return {"bundle_id": bundle_id, "bundle_manifest": str(bundle_manifest_path), "n_profiles": len(consensus_profiles)}
