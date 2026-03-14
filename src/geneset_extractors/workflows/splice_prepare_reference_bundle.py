from __future__ import annotations

import csv
import gzip
import json
import math
from pathlib import Path
from statistics import mean

from geneset_extractors.hashing import sha256_file
from geneset_extractors.extractors.splicing.splice_event_diff_workflow import _normalize_event_type



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

    organism = str(args.organism).strip().lower()
    alias_records: dict[str, dict[str, object]] = {}
    ubiquity_stats: dict[str, dict[str, object]] = {}
    impact_records: dict[str, dict[str, object]] = {}
    source_file_records: list[dict[str, object]] = []
    total_samples_seen: set[str] = set()

    for source_row in source_rows:
        source_file = Path(_clean(source_row.get("path")))
        if not source_file.is_absolute() and not source_file.exists():
            source_file = (sources_path.parent / source_file).resolve()
        source_dataset = _clean(source_row.get("source_dataset")) or source_file.stem
        with source_file.open("r", encoding="utf-8", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        present_by_event: dict[str, set[str]] = {}
        for row in rows:
            input_event_key = _clean(row.get("input_event_key"))
            canonical = _clean(row.get("canonical_event_key"))
            if not input_event_key or not canonical:
                continue
            gene_id = _clean(row.get("gene_id"))
            gene_symbol = _clean(row.get("gene_symbol")) or gene_id
            event_type = _normalize_event_type(_clean(row.get("event_type")))
            sample_id = _clean(row.get("sample_id"))
            psi = _parse_float(row.get("psi"))
            read_support = _parse_float(row.get("read_support"))
            annotation_status = _clean(row.get("annotation_status"))
            total_samples_seen.add(sample_id)
            alias_records[input_event_key] = {
                "input_event_key": input_event_key,
                "canonical_event_key": canonical,
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "event_type": event_type,
                "chrom": _clean(row.get("chrom")),
                "start": _clean(row.get("start")),
                "end": _clean(row.get("end")),
                "strand": _clean(row.get("strand")),
                "source_dataset": source_dataset,
            }
            if psi is not None and (read_support is None or read_support >= float(args.min_ref_read_support)):
                present_by_event.setdefault(canonical, set()).add(sample_id)
            stats = ubiquity_stats.setdefault(canonical, {
                "canonical_event_key": canonical,
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "event_type": event_type,
                "datasets": set(),
                "df_ref": 0,
            })
            stats["datasets"].add(source_dataset)
            impact_records.setdefault(canonical, {
                "canonical_event_key": canonical,
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "event_type": event_type,
                "annotation_status": annotation_status,
            })
        for canonical, sample_ids in present_by_event.items():
            ubiquity_stats[canonical]["df_ref"] = int(ubiquity_stats[canonical]["df_ref"]) + len(sample_ids)
        source_file_records.append({
            "path": str(source_file),
            "sha256": sha256_file(source_file),
            "source_dataset": source_dataset,
            "n_rows": len(rows),
            "n_unique_events": len(present_by_event),
        })

    n_samples_ref = max(1, len(total_samples_seen))
    alias_rows = [alias_records[key] for key in sorted(alias_records)]
    ubiquity_rows: list[dict[str, object]] = []
    idf_values: list[float] = []
    for canonical in sorted(ubiquity_stats):
        stats = ubiquity_stats[canonical]
        df_ref = int(stats["df_ref"])
        idf_ref = math.log1p((n_samples_ref + 1.0) / (df_ref + 1.0))
        idf_values.append(idf_ref)
        ubiquity_rows.append({
            "canonical_event_key": canonical,
            "gene_id": stats["gene_id"],
            "gene_symbol": stats["gene_symbol"],
            "event_type": stats["event_type"],
            "n_samples_ref": n_samples_ref,
            "df_ref": df_ref,
            "fraction_ref": float(df_ref) / float(n_samples_ref),
            "idf_ref": idf_ref,
            "n_datasets_ref": len(stats["datasets"]),
        })

    impact_rows: list[dict[str, object]] = []
    for canonical in sorted(impact_records):
        record = impact_records[canonical]
        raw, evidence = _impact_defaults(str(record["event_type"]), str(record["annotation_status"]))
        impact_rows.append({
            "canonical_event_key": canonical,
            "gene_id": record["gene_id"],
            "gene_symbol": record["gene_symbol"],
            "event_type": record["event_type"],
            "impact_weight_raw": raw,
            "impact_evidence": evidence,
            "annotation_status": record["annotation_status"],
        })

    alias_filename = f"splice_event_aliases_{organism}_v1.tsv.gz"
    ubiquity_filename = f"splice_event_ubiquity_{organism}_v1.tsv.gz"
    impact_filename = f"splice_event_impact_{organism}_v1.tsv.gz"
    alias_path = out_dir / alias_filename
    ubiquity_path = out_dir / ubiquity_filename
    impact_path = out_dir / impact_filename
    provenance_path = out_dir / "bundle_provenance.json"
    local_manifest_path = out_dir / "local_resources_manifest.json"

    _write_tsv_gz(alias_path, ["input_event_key", "canonical_event_key", "gene_id", "gene_symbol", "event_type", "chrom", "start", "end", "strand", "source_dataset"], alias_rows)
    _write_tsv_gz(ubiquity_path, ["canonical_event_key", "gene_id", "gene_symbol", "event_type", "n_samples_ref", "df_ref", "fraction_ref", "idf_ref", "n_datasets_ref"], ubiquity_rows)
    _write_tsv_gz(impact_path, ["canonical_event_key", "gene_id", "gene_symbol", "event_type", "impact_weight_raw", "impact_evidence", "annotation_status"], impact_rows)

    provenance = {
        "bundle_id": args.bundle_id,
        "organism": organism,
        "n_sources": len(source_rows),
        "n_alias_rows": len(alias_rows),
        "n_canonical_events": len(ubiquity_rows),
        "idf_ref_summary": {"min": min(idf_values) if idf_values else 0.0, "mean": mean(idf_values) if idf_values else 0.0, "max": max(idf_values) if idf_values else 0.0},
        "source_files": source_file_records,
        "outputs": {"alias_table": alias_filename, "ubiquity_table": ubiquity_filename, "impact_table": impact_filename},
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
            }
        ],
        "presets": {
            f"splice_event_bundle_{organism}_v1": [
                f"splice_event_aliases_{organism}_v1",
                f"splice_event_ubiquity_{organism}_v1",
                f"splice_event_impact_{organism}_v1"
            ]
        }
    }
    local_manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "workflow": "splice_prepare_reference_bundle",
        "bundle_id": args.bundle_id,
        "out_dir": str(out_dir),
        "n_canonical_events": len(ubiquity_rows),
        "n_alias_rows": len(alias_rows),
        "n_impact_rows": len(impact_rows),
    }
