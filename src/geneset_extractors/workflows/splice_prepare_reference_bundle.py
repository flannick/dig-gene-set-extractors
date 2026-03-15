from __future__ import annotations

import csv
import gzip
import json
import math
from pathlib import Path
from statistics import mean, median

from geneset_extractors.extractors.splicing.splice_event_diff_workflow import _normalize_event_type
from geneset_extractors.hashing import sha256_file


def _clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _parse_float(value: object) -> float | None:
    text = _clean(value)
    if not text or text.lower() in {"na", "nan", "none", "null"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _parse_csv_set(value: object) -> set[str]:
    text = _clean(value)
    if not text:
        return set()
    return {token.strip() for token in text.split(",") if token.strip()}


def _write_tsv_gz(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _impact_defaults(event_type: str, annotation_status: str) -> tuple[float, float]:
    raw = 1.0
    evidence = 0.2
    et = _normalize_event_type(event_type)
    ann = annotation_status.lower()
    if et in {"exon_skip", "alt_donor", "alt_acceptor", "mutually_exclusive_exon"}:
        raw = 1.12
        evidence = 0.45
    elif et == "retained_intron":
        raw = 0.9
        evidence = 0.4
    if "coding" in ann:
        raw = max(raw, 1.15)
        evidence = max(evidence, 0.55)
    if "novel" in ann or "unannotated" in ann:
        raw = min(raw, 0.95)
        evidence = max(evidence, 0.3)
    return raw, evidence


def _event_key_namespace(canonicalization_status: str, canonicalization_confidence: str, value: str) -> str:
    explicit = _clean(value)
    if explicit:
        return explicit
    status = _clean(canonicalization_status)
    confidence = _clean(canonicalization_confidence).lower()
    if status == "coordinate_canonical" or confidence == "high":
        return "global_coordinate"
    if status == "source_family_stable_id" or confidence == "medium":
        return "tcga_spliceseq_asid"
    return "raw_id_fallback"


def _bundle_event_key(canonical_event_key: str, source_dataset: str, canonicalization_confidence: str) -> str:
    if canonicalization_confidence in {"high", "medium"}:
        return canonical_event_key
    return f"lowconf::{source_dataset}::{canonical_event_key}"


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sources_path = Path(args.sources_tsv)
    with sources_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        source_rows = list(reader)
    if not source_rows:
        raise ValueError("sources_tsv must contain at least one row")
    if "path" not in reader.fieldnames:
        raise ValueError("sources_tsv must contain a path column")
    excluded_datasets = _parse_csv_set(getattr(args, "exclude_source_datasets", ""))
    if excluded_datasets:
        source_rows = [
            row
            for row in source_rows
            if (_clean(row.get("source_dataset")) or Path(_clean(row.get("path"))).stem) not in excluded_datasets
        ]
        if not source_rows:
            raise ValueError("All source rows were excluded by --exclude_source_datasets.")

    organism = str(args.organism).strip().lower()
    alias_records: dict[tuple[str, str], dict[str, object]] = {}
    alias_conflicts: dict[str, int] = {}
    ubiquity_stats: dict[str, dict[str, object]] = {}
    ubiquity_by_dataset_stats: dict[tuple[str, str], dict[str, object]] = {}
    impact_records: dict[str, dict[str, object]] = {}
    gene_burden: dict[str, dict[str, object]] = {}
    gene_burden_by_dataset: dict[tuple[str, str], dict[str, object]] = {}
    source_file_records: list[dict[str, object]] = []
    total_samples_seen: set[str] = set()
    total_source_datasets: set[str] = set()
    canonicalization_status_counts: dict[str, int] = {
        "coordinate_canonical": 0,
        "source_family_stable_id": 0,
        "raw_id_fallback": 0,
    }
    canonicalization_confidence_counts: dict[str, int] = {"high": 0, "medium": 0, "low": 0}

    for source_row in source_rows:
        source_file = Path(_clean(source_row.get("path")))
        if not source_file.is_absolute() and not source_file.exists():
            source_file = (sources_path.parent / source_file).resolve()
        source_dataset = _clean(source_row.get("source_dataset")) or source_file.stem
        total_source_datasets.add(source_dataset)
        with source_file.open("r", encoding="utf-8", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        present_by_event: dict[str, set[str]] = {}
        unique_events_by_gene: dict[str, set[str]] = {}
        unique_groups_by_gene: dict[str, set[str]] = {}
        unique_groups_by_gene_study: dict[str, dict[str, set[str]]] = {}
        studies_by_gene: dict[str, set[str]] = {}
        high_conf_studies_by_gene: dict[str, set[str]] = {}
        unique_high_events_by_gene: dict[str, set[str]] = {}
        unique_medium_events_by_gene: dict[str, set[str]] = {}
        unique_low_events_by_gene: dict[str, set[str]] = {}
        for row in rows:
            input_event_key = _clean(row.get("input_event_key"))
            canonical = _clean(row.get("canonical_event_key"))
            if not input_event_key or not canonical:
                continue
            canonicalization_status = _clean(row.get("canonicalization_status")) or "raw_id_fallback"
            canonicalization_confidence = _clean(row.get("canonicalization_confidence")).lower() or ("high" if canonicalization_status == "coordinate_canonical" else "low")
            event_key_namespace = _event_key_namespace(
                canonicalization_status,
                canonicalization_confidence,
                _clean(row.get("event_key_namespace")),
            )
            canonicalization_status_counts[canonicalization_status] = canonicalization_status_counts.get(canonicalization_status, 0) + 1
            canonicalization_confidence_counts[canonicalization_confidence] = canonicalization_confidence_counts.get(canonicalization_confidence, 0) + 1
            bundle_event_key = _bundle_event_key(canonical, source_dataset, canonicalization_confidence)
            gene_id = _clean(row.get("gene_id"))
            gene_symbol = _clean(row.get("gene_symbol")) or gene_id
            event_type = _normalize_event_type(_clean(row.get("event_type")))
            event_group = _clean(row.get("event_group")) or canonical
            sample_id = _clean(row.get("sample_id"))
            psi = _parse_float(row.get("psi"))
            read_support = _parse_float(row.get("read_support"))
            annotation_status = _clean(row.get("annotation_status"))
            total_samples_seen.add(sample_id)

            alias_key = (input_event_key, canonicalization_confidence)
            alias_record = {
                "input_event_key": input_event_key,
                "canonical_event_key": bundle_event_key,
                "canonicalization_status": canonicalization_status,
                "canonicalization_confidence": canonicalization_confidence,
                "event_key_namespace": event_key_namespace,
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "event_type": event_type,
                "chrom": _clean(row.get("chrom")),
                "start": _clean(row.get("start")),
                "end": _clean(row.get("end")),
                "strand": _clean(row.get("strand")),
                "source_dataset": source_dataset,
            }
            if alias_key in alias_records:
                prev = alias_records[alias_key]
                same_mapping = (
                    str(prev.get("canonical_event_key")) == alias_record["canonical_event_key"]
                    and str(prev.get("gene_symbol")) == alias_record["gene_symbol"]
                )
                if not same_mapping and canonicalization_confidence in {"low", "medium"}:
                    alias_conflicts[input_event_key] = alias_conflicts.get(input_event_key, 0) + 1
                    del alias_records[alias_key]
                elif same_mapping:
                    alias_records[alias_key] = prev
            elif alias_conflicts.get(input_event_key, 0) == 0:
                alias_records[alias_key] = alias_record

            if psi is not None and (read_support is None or read_support >= float(args.min_ref_read_support)):
                present_by_event.setdefault(bundle_event_key, set()).add(sample_id)
            dataset_stats = ubiquity_by_dataset_stats.setdefault(
                (bundle_event_key, source_dataset),
                {
                    "canonical_event_key": bundle_event_key,
                    "source_dataset": source_dataset,
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "event_type": event_type,
                    "canonicalization_status": canonicalization_status,
                    "canonicalization_confidence": canonicalization_confidence,
                    "event_key_namespace": event_key_namespace,
                    "df_ref": 0,
                    "samples": set(),
                },
            )
            if sample_id:
                dataset_stats["samples"].add(sample_id)
            stats = ubiquity_stats.setdefault(
                bundle_event_key,
                {
                    "canonical_event_key": bundle_event_key,
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "event_type": event_type,
                    "canonicalization_status": canonicalization_status,
                    "canonicalization_confidence": canonicalization_confidence,
                    "event_key_namespace": event_key_namespace,
                    "datasets": set(),
                    "samples": set(),
                    "df_ref": 0,
                },
            )
            stats["datasets"].add(source_dataset)
            if sample_id:
                stats["samples"].add(sample_id)
            impact_records.setdefault(
                bundle_event_key,
                {
                    "canonical_event_key": bundle_event_key,
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "event_type": event_type,
                    "annotation_status": annotation_status,
                    "canonicalization_status": canonicalization_status,
                    "canonicalization_confidence": canonicalization_confidence,
                    "event_key_namespace": event_key_namespace,
                    "datasets": set(),
                },
            )
            impact_records[bundle_event_key]["datasets"].add(source_dataset)
            if gene_symbol:
                unique_events_by_gene.setdefault(gene_symbol, set()).add(bundle_event_key)
                unique_groups_by_gene.setdefault(gene_symbol, set()).add(event_group or bundle_event_key)
                unique_groups_by_gene_study.setdefault(gene_symbol, {}).setdefault(source_dataset, set()).add(event_group or bundle_event_key)
                studies_by_gene.setdefault(gene_symbol, set()).add(source_dataset)
                if canonicalization_confidence == "high":
                    unique_high_events_by_gene.setdefault(gene_symbol, set()).add(bundle_event_key)
                    high_conf_studies_by_gene.setdefault(gene_symbol, set()).add(source_dataset)
                elif canonicalization_confidence == "medium":
                    unique_medium_events_by_gene.setdefault(gene_symbol, set()).add(bundle_event_key)
                else:
                    unique_low_events_by_gene.setdefault(gene_symbol, set()).add(bundle_event_key)
                burden_by_dataset = gene_burden_by_dataset.setdefault(
                    (gene_symbol, source_dataset),
                    {
                        "gene_symbol": gene_symbol,
                        "source_dataset": source_dataset,
                        "n_canonical_events_ref": set(),
                        "n_high_confidence_events_ref": set(),
                        "n_medium_confidence_events_ref": set(),
                        "n_low_confidence_events_ref": set(),
                        "n_unique_event_groups_ref": set(),
                    },
                )
                burden_by_dataset["n_canonical_events_ref"].add(bundle_event_key)
                burden_by_dataset["n_unique_event_groups_ref"].add(event_group or bundle_event_key)
                if canonicalization_confidence == "high":
                    burden_by_dataset["n_high_confidence_events_ref"].add(bundle_event_key)
                elif canonicalization_confidence == "medium":
                    burden_by_dataset["n_medium_confidence_events_ref"].add(bundle_event_key)
                else:
                    burden_by_dataset["n_low_confidence_events_ref"].add(bundle_event_key)
        for canonical_key, sample_ids in present_by_event.items():
            ubiquity_stats[canonical_key]["df_ref"] = int(ubiquity_stats[canonical_key]["df_ref"]) + len(sample_ids)
            if (canonical_key, source_dataset) in ubiquity_by_dataset_stats:
                ubiquity_by_dataset_stats[(canonical_key, source_dataset)]["df_ref"] = len(sample_ids)
        for gene_symbol, bundle_keys in unique_events_by_gene.items():
            burden = gene_burden.setdefault(
                gene_symbol,
                {
                    "gene_symbol": gene_symbol,
                    "n_canonical_events_ref": set(),
                    "n_high_confidence_events_ref": set(),
                    "n_medium_confidence_events_ref": set(),
                    "n_low_confidence_events_ref": set(),
                    "n_unique_event_groups_ref": set(),
                    "studies": set(),
                    "high_conf_studies": set(),
                    "study_unique_groups": {},
                },
            )
            burden["n_canonical_events_ref"].update(bundle_keys)
            burden["n_high_confidence_events_ref"].update(unique_high_events_by_gene.get(gene_symbol, set()))
            burden["n_medium_confidence_events_ref"].update(unique_medium_events_by_gene.get(gene_symbol, set()))
            burden["n_low_confidence_events_ref"].update(unique_low_events_by_gene.get(gene_symbol, set()))
            burden["n_unique_event_groups_ref"].update(unique_groups_by_gene.get(gene_symbol, set()))
            burden["studies"].update(studies_by_gene.get(gene_symbol, set()))
            burden["high_conf_studies"].update(high_conf_studies_by_gene.get(gene_symbol, set()))
            for study_id, group_set in unique_groups_by_gene_study.get(gene_symbol, {}).items():
                burden["study_unique_groups"].setdefault(study_id, set()).update(group_set)
        source_file_records.append(
            {
                "path": str(source_file),
                "sha256": sha256_file(source_file),
                "source_dataset": source_dataset,
                "n_rows": len(rows),
                "n_unique_events": len(present_by_event),
            }
        )

    n_samples_ref = max(1, len(total_samples_seen))
    n_datasets_total = max(1, len(total_source_datasets))
    alias_rows = [alias_records[key] for key in sorted(alias_records)]
    ubiquity_rows: list[dict[str, object]] = []
    ubiquity_by_dataset_rows: list[dict[str, object]] = []
    idf_values: list[float] = []
    for canonical in sorted(ubiquity_stats):
        stats = ubiquity_stats[canonical]
        df_ref = int(stats["df_ref"])
        n_datasets_ref = len(stats["datasets"])
        idf_ref = math.log1p((n_datasets_total + 1.0) / (n_datasets_ref + 1.0))
        idf_values.append(idf_ref)
        ubiquity_rows.append(
            {
                "canonical_event_key": canonical,
                "gene_id": stats["gene_id"],
                "gene_symbol": stats["gene_symbol"],
                "event_type": stats["event_type"],
                "canonicalization_status": stats["canonicalization_status"],
                "canonicalization_confidence": stats["canonicalization_confidence"],
                "event_key_namespace": stats["event_key_namespace"],
                "n_samples_ref": n_samples_ref,
                "df_ref": df_ref,
                "fraction_ref": float(df_ref) / float(n_samples_ref),
                "idf_ref": idf_ref,
                "n_datasets_ref": n_datasets_ref,
                "fraction_datasets_ref": float(n_datasets_ref) / float(n_datasets_total),
            }
        )
    for (_canonical, _dataset) in sorted(ubiquity_by_dataset_stats):
        stats = ubiquity_by_dataset_stats[(_canonical, _dataset)]
        ubiquity_by_dataset_rows.append(
            {
                "canonical_event_key": stats["canonical_event_key"],
                "source_dataset": stats["source_dataset"],
                "gene_id": stats["gene_id"],
                "gene_symbol": stats["gene_symbol"],
                "event_type": stats["event_type"],
                "canonicalization_status": stats["canonicalization_status"],
                "canonicalization_confidence": stats["canonicalization_confidence"],
                "event_key_namespace": stats["event_key_namespace"],
                "n_samples_ref": len(stats["samples"]),
                "df_ref": int(stats["df_ref"]),
            }
        )

    impact_rows: list[dict[str, object]] = []
    for canonical in sorted(impact_records):
        record = impact_records[canonical]
        raw, evidence = _impact_defaults(str(record["event_type"]), str(record["annotation_status"]))
        impact_rows.append(
            {
                "canonical_event_key": canonical,
                "gene_id": record["gene_id"],
                "gene_symbol": record["gene_symbol"],
                "event_type": record["event_type"],
                "canonicalization_status": record["canonicalization_status"],
                "canonicalization_confidence": record["canonicalization_confidence"],
                "event_key_namespace": record["event_key_namespace"],
                "impact_weight_raw": raw,
                "impact_evidence": evidence,
                "annotation_status": record["annotation_status"],
                "n_datasets_ref": len(record["datasets"]),
            }
        )

    burden_rows: list[dict[str, object]] = []
    burden_by_dataset_rows: list[dict[str, object]] = []
    for gene_symbol in sorted(gene_burden):
        record = gene_burden[gene_symbol]
        study_group_counts = [
            len(group_set)
            for group_set in record["study_unique_groups"].values()
            if group_set
        ]
        total_groups = len(record["n_unique_event_groups_ref"])
        low_fraction = (
            float(len(record["n_low_confidence_events_ref"])) / float(max(1, len(record["n_canonical_events_ref"])))
        )
        burden_rows.append(
            {
                "gene_symbol": gene_symbol,
                "n_canonical_events_ref": len(record["n_canonical_events_ref"]),
                "n_high_confidence_events_ref": len(record["n_high_confidence_events_ref"]),
                "n_medium_confidence_events_ref": len(record["n_medium_confidence_events_ref"]),
                "n_low_confidence_events_ref": len(record["n_low_confidence_events_ref"]),
                "n_unique_event_groups_ref": total_groups,
                "n_studies_ref": len(record["studies"]),
                "n_studies_high_confidence_ref": len(record["high_conf_studies"]),
                "fraction_low_confidence_events_ref": low_fraction,
                "median_unique_groups_per_study": float(median(study_group_counts)) if study_group_counts else 0.0,
            }
        )
    for (_gene_symbol, _dataset) in sorted(gene_burden_by_dataset):
        record = gene_burden_by_dataset[(_gene_symbol, _dataset)]
        burden_by_dataset_rows.append(
            {
                "gene_symbol": record["gene_symbol"],
                "source_dataset": record["source_dataset"],
                "n_canonical_events_ref": len(record["n_canonical_events_ref"]),
                "n_high_confidence_events_ref": len(record["n_high_confidence_events_ref"]),
                "n_medium_confidence_events_ref": len(record["n_medium_confidence_events_ref"]),
                "n_low_confidence_events_ref": len(record["n_low_confidence_events_ref"]),
                "n_unique_event_groups_ref": len(record["n_unique_event_groups_ref"]),
            }
        )

    alias_filename = f"splice_event_aliases_{organism}_v1.tsv.gz"
    ubiquity_filename = f"splice_event_ubiquity_{organism}_v1.tsv.gz"
    ubiquity_by_dataset_filename = f"splice_event_ubiquity_by_dataset_{organism}_v1.tsv.gz"
    impact_filename = f"splice_event_impact_{organism}_v1.tsv.gz"
    burden_filename = f"splice_gene_event_burden_{organism}_v1.tsv.gz"
    burden_by_dataset_filename = f"splice_gene_event_burden_by_dataset_{organism}_v1.tsv.gz"
    alias_path = out_dir / alias_filename
    ubiquity_path = out_dir / ubiquity_filename
    ubiquity_by_dataset_path = out_dir / ubiquity_by_dataset_filename
    impact_path = out_dir / impact_filename
    burden_path = out_dir / burden_filename
    burden_by_dataset_path = out_dir / burden_by_dataset_filename
    provenance_path = out_dir / "bundle_provenance.json"
    local_manifest_path = out_dir / "local_resources_manifest.json"

    _write_tsv_gz(
        alias_path,
        [
            "input_event_key",
            "canonical_event_key",
            "canonicalization_status",
            "canonicalization_confidence",
            "event_key_namespace",
            "gene_id",
            "gene_symbol",
            "event_type",
            "chrom",
            "start",
            "end",
            "strand",
            "source_dataset",
        ],
        alias_rows,
    )
    _write_tsv_gz(
        ubiquity_path,
        [
            "canonical_event_key",
            "gene_id",
            "gene_symbol",
            "event_type",
            "canonicalization_status",
            "canonicalization_confidence",
            "event_key_namespace",
            "n_samples_ref",
            "df_ref",
            "fraction_ref",
            "idf_ref",
            "n_datasets_ref",
            "fraction_datasets_ref",
        ],
        ubiquity_rows,
    )
    _write_tsv_gz(
        ubiquity_by_dataset_path,
        [
            "canonical_event_key",
            "source_dataset",
            "gene_id",
            "gene_symbol",
            "event_type",
            "canonicalization_status",
            "canonicalization_confidence",
            "event_key_namespace",
            "n_samples_ref",
            "df_ref",
        ],
        ubiquity_by_dataset_rows,
    )
    _write_tsv_gz(
        impact_path,
        [
            "canonical_event_key",
            "gene_id",
            "gene_symbol",
            "event_type",
            "canonicalization_status",
            "canonicalization_confidence",
            "event_key_namespace",
            "impact_weight_raw",
            "impact_evidence",
            "annotation_status",
            "n_datasets_ref",
        ],
        impact_rows,
    )
    _write_tsv_gz(
        burden_path,
        [
            "gene_symbol",
            "n_canonical_events_ref",
            "n_high_confidence_events_ref",
            "n_medium_confidence_events_ref",
            "n_low_confidence_events_ref",
            "n_unique_event_groups_ref",
            "n_studies_ref",
            "n_studies_high_confidence_ref",
            "fraction_low_confidence_events_ref",
            "median_unique_groups_per_study",
        ],
        burden_rows,
    )
    _write_tsv_gz(
        burden_by_dataset_path,
        [
            "gene_symbol",
            "source_dataset",
            "n_canonical_events_ref",
            "n_high_confidence_events_ref",
            "n_medium_confidence_events_ref",
            "n_low_confidence_events_ref",
            "n_unique_event_groups_ref",
        ],
        burden_by_dataset_rows,
    )

    provenance = {
        "bundle_id": args.bundle_id,
        "organism": organism,
        "n_sources": len(source_rows),
        "n_source_datasets": len(total_source_datasets),
        "excluded_source_datasets": sorted(excluded_datasets),
        "n_alias_rows": len(alias_rows),
        "n_canonical_events": len(ubiquity_rows),
        "n_genes_with_burden": len(burden_rows),
        "n_event_ubiquity_by_dataset_rows": len(ubiquity_by_dataset_rows),
        "n_gene_burden_by_dataset_rows": len(burden_by_dataset_rows),
        "canonicalization_status_counts": canonicalization_status_counts,
        "canonicalization_confidence_counts": canonicalization_confidence_counts,
        "fraction_high_confidence_events": (
            float(canonicalization_confidence_counts.get("high", 0)) / float(max(1, sum(canonicalization_confidence_counts.values())))
        ),
        "fraction_medium_confidence_events": (
            float(canonicalization_confidence_counts.get("medium", 0)) / float(max(1, sum(canonicalization_confidence_counts.values())))
        ),
        "fraction_low_confidence_events": (
            float(canonicalization_confidence_counts.get("low", 0)) / float(max(1, sum(canonicalization_confidence_counts.values())))
        ),
        "low_confidence_alias_conflicts_dropped": sum(alias_conflicts.values()),
        "idf_ref_summary": {
            "min": min(idf_values) if idf_values else 0.0,
            "mean": mean(idf_values) if idf_values else 0.0,
            "max": max(idf_values) if idf_values else 0.0,
        },
        "idf_ref_basis": "study_level_dataset_frequency",
        "source_files": source_file_records,
        "outputs": {
            "alias_table": alias_filename,
            "ubiquity_table": ubiquity_filename,
            "ubiquity_by_dataset_table": ubiquity_by_dataset_filename,
            "impact_table": impact_filename,
            "gene_burden_table": burden_filename,
            "gene_burden_by_dataset_table": burden_by_dataset_filename,
        },
    }
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    manifest = {
        "resources": [
            {
                "id": f"splice_event_aliases_{organism}_v1",
                "description": f"Canonical splice-event alias table for {organism}",
                "provider": "local_bundle",
                "stable_id": f"local:splice_event_aliases_{organism}_v1",
                "version": "v1",
                "genome_build": organism,
                "filename": alias_filename,
                "url": "",
                "sha256": sha256_file(alias_path),
                "license": "Local bundle resource; see bundle_provenance.json",
            },
            {
                "id": f"splice_event_ubiquity_{organism}_v1",
                "description": f"Splice-event ubiquity prior for {organism}",
                "provider": "local_bundle",
                "stable_id": f"local:splice_event_ubiquity_{organism}_v1",
                "version": "v1",
                "genome_build": organism,
                "filename": ubiquity_filename,
                "url": "",
                "sha256": sha256_file(ubiquity_path),
                "license": "Local bundle resource; see bundle_provenance.json",
            },
            {
                "id": f"splice_event_ubiquity_by_dataset_{organism}_v1",
                "description": f"Dataset-stratified splice-event ubiquity summary for {organism}",
                "provider": "local_bundle",
                "stable_id": f"local:splice_event_ubiquity_by_dataset_{organism}_v1",
                "version": "v1",
                "genome_build": organism,
                "filename": ubiquity_by_dataset_filename,
                "url": "",
                "sha256": sha256_file(ubiquity_by_dataset_path),
                "license": "Local bundle resource; see bundle_provenance.json",
            },
            {
                "id": f"splice_event_impact_{organism}_v1",
                "description": f"Conservative splice-event impact prior for {organism}",
                "provider": "local_bundle",
                "stable_id": f"local:splice_event_impact_{organism}_v1",
                "version": "v1",
                "genome_build": organism,
                "filename": impact_filename,
                "url": "",
                "sha256": sha256_file(impact_path),
                "license": "Local bundle resource; see bundle_provenance.json",
            },
            {
                "id": f"splice_gene_event_burden_{organism}_v1",
                "description": f"Per-gene splicing event burden reference counts for {organism}",
                "provider": "local_bundle",
                "stable_id": f"local:splice_gene_event_burden_{organism}_v1",
                "version": "v1",
                "genome_build": organism,
                "filename": burden_filename,
                "url": "",
                "sha256": sha256_file(burden_path),
                "license": "Local bundle resource; see bundle_provenance.json",
            },
            {
                "id": f"splice_gene_event_burden_by_dataset_{organism}_v1",
                "description": f"Dataset-stratified per-gene splicing burden counts for {organism}",
                "provider": "local_bundle",
                "stable_id": f"local:splice_gene_event_burden_by_dataset_{organism}_v1",
                "version": "v1",
                "genome_build": organism,
                "filename": burden_by_dataset_filename,
                "url": "",
                "sha256": sha256_file(burden_by_dataset_path),
                "license": "Local bundle resource; see bundle_provenance.json",
            },
        ],
        "presets": {
            f"splice_event_bundle_{organism}_v1": [
                f"splice_event_aliases_{organism}_v1",
                f"splice_event_ubiquity_{organism}_v1",
                f"splice_event_ubiquity_by_dataset_{organism}_v1",
                f"splice_event_impact_{organism}_v1",
                f"splice_gene_event_burden_{organism}_v1",
                f"splice_gene_event_burden_by_dataset_{organism}_v1",
            ]
        },
    }
    local_manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "workflow": "splice_prepare_reference_bundle",
        "bundle_id": args.bundle_id,
        "out_dir": str(out_dir),
        "n_canonical_events": len(ubiquity_rows),
        "n_alias_rows": len(alias_rows),
        "n_impact_rows": len(impact_rows),
        "n_gene_burden_rows": len(burden_rows),
        "n_event_ubiquity_by_dataset_rows": len(ubiquity_by_dataset_rows),
        "n_gene_burden_by_dataset_rows": len(burden_by_dataset_rows),
        "n_source_datasets": len(total_source_datasets),
    }
