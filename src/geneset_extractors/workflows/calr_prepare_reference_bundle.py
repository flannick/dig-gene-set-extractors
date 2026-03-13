from __future__ import annotations

import csv
import gzip
import hashlib
import json
from contextlib import ExitStack
from importlib.resources import as_file
from pathlib import Path
from statistics import mean
import tarfile

from geneset_extractors.extractors.calorimetry.ontology import packaged_resource_path
from geneset_extractors.extractors.calorimetry.orthology import packaged_resource_path as packaged_orthology_resource_path
from geneset_extractors.extractors.morphology.io import write_bundle_manifest


DEFAULT_TERM_TEMPLATES = "calr_term_templates_mouse_v1.tsv"
DEFAULT_GENE_EDGES = "calr_phenotype_gene_edges_mouse_v1.tsv"
DEFAULT_TERM_HIERARCHY = "calr_term_hierarchy_mouse_v1.tsv"


def _open_text(path: str | Path):
    p = Path(path)
    if p.suffix.lower() == ".gz":
        return gzip.open(p, "rt", encoding="utf-8")
    return p.open("r", encoding="utf-8")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _read_rows(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing header in {path}")
        fieldnames = [str(name) for name in reader.fieldnames]
        rows = [{str(k): str(v or "") for k, v in row.items()} for row in reader]
    return fieldnames, rows


def _write_tsv_gz(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _copy_to_tsv_gz(src: str | Path, dest: Path) -> tuple[int, str]:
    fieldnames, rows = _read_rows(src)
    _write_tsv_gz(dest, fieldnames, rows)
    return len(rows), str(src)


def _derive_feature_schema_and_stats(reference_profiles_tsv: str | Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    fieldnames, rows = _read_rows(reference_profiles_tsv)
    if len(fieldnames) < 2:
        raise ValueError("reference_profiles_tsv must have an id column plus at least one feature column")
    feature_names = fieldnames[1:]
    schema_rows = [{"feature_name": feature_name} for feature_name in feature_names]
    stats_rows: list[dict[str, object]] = []
    for feature_name in feature_names:
        values = []
        for row in rows:
            text = str(row.get(feature_name, "")).strip()
            if not text or text.lower() in {"na", "nan", "null", "none"}:
                continue
            try:
                values.append(float(text))
            except ValueError:
                continue
        center = float(mean(values)) if values else 0.0
        scale = float(max(1e-8, (sum((value - center) ** 2 for value in values) / float(len(values) or 1)) ** 0.5))
        stats_rows.append({"feature_name": feature_name, "center": center, "scale": scale})
    return schema_rows, stats_rows


def _write_distribution_artifact(*, dist_dir: Path, bundle_id: str, included_files: list[Path]) -> dict[str, str]:
    dist_dir.mkdir(parents=True, exist_ok=True)
    tarball_path = dist_dir / f"dig-gene-set-extractors-{bundle_id}.tar.gz"
    with tarfile.open(tarball_path, "w:gz") as tf:
        for path in included_files:
            tf.add(path, arcname=f"bundle/{path.name}")
    tar_sha = _sha256_file(tarball_path)
    sha_path = dist_dir / "SHA256SUMS.txt"
    sha_path.write_text(f"{tar_sha}  {tarball_path.name}\n", encoding="utf-8")
    manifest_path = dist_dir / "distribution_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "bundle_id": bundle_id,
                "tarball": tarball_path.name,
                "sha256sums": sha_path.name,
                "bundle_root_inside_tarball": "bundle",
                "included_files": [path.name for path in included_files],
                "tarball_sha256": tar_sha,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "distribution_dir": str(dist_dir),
        "tarball": str(tarball_path),
        "sha256sums": str(sha_path),
        "distribution_manifest": str(manifest_path),
    }


def run(args) -> dict[str, object]:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bundle_id = str(args.bundle_id or "").strip() or "calorimetry_reference_bundle_v1"

    with ExitStack() as stack:
        term_templates_src = args.term_templates_tsv
        phenotype_gene_edges_src = args.phenotype_gene_edges_tsv
        term_hierarchy_src = args.term_hierarchy_tsv
        ortholog_src = getattr(args, "mouse_human_orthologs_tsv", None)
        if not term_templates_src:
            term_templates_src = str(stack.enter_context(as_file(packaged_resource_path(DEFAULT_TERM_TEMPLATES))))
        if not phenotype_gene_edges_src:
            phenotype_gene_edges_src = str(stack.enter_context(as_file(packaged_resource_path(DEFAULT_GENE_EDGES))))
        if not term_hierarchy_src and bool(args.include_packaged_term_hierarchy):
            term_hierarchy_src = str(stack.enter_context(as_file(packaged_resource_path(DEFAULT_TERM_HIERARCHY))))
        if str(args.organism).strip().lower() == "mouse" and not ortholog_src:
            ortholog_src = str(stack.enter_context(as_file(packaged_orthology_resource_path("calr_mouse_human_orthologs_v1.tsv.gz"))))

        profiles_out = out_dir / "reference_profiles.tsv.gz"
        metadata_out = out_dir / "reference_metadata.tsv.gz"
        schema_out = out_dir / "feature_schema.tsv.gz"
        stats_out = out_dir / "feature_stats.tsv.gz"
        templates_out = out_dir / "term_templates.tsv.gz"
        edges_out = out_dir / "phenotype_gene_edges.tsv.gz"
        hierarchy_out = out_dir / "term_hierarchy.tsv.gz"
        ortholog_out = out_dir / "mouse_human_orthologs.tsv.gz"

        n_profiles, _ = _copy_to_tsv_gz(args.reference_profiles_tsv, profiles_out)
        n_metadata, _ = _copy_to_tsv_gz(args.reference_metadata_tsv, metadata_out)
        n_template_rows, _ = _copy_to_tsv_gz(term_templates_src, templates_out)
        n_edge_rows, _ = _copy_to_tsv_gz(phenotype_gene_edges_src, edges_out)
        n_hierarchy_rows = 0
        if term_hierarchy_src:
            n_hierarchy_rows, _ = _copy_to_tsv_gz(term_hierarchy_src, hierarchy_out)
        n_ortholog_rows = 0
        if ortholog_src:
            n_ortholog_rows, _ = _copy_to_tsv_gz(ortholog_src, ortholog_out)

        if args.feature_schema_tsv:
            n_schema_rows, _ = _copy_to_tsv_gz(args.feature_schema_tsv, schema_out)
        else:
            schema_rows, stats_rows = _derive_feature_schema_and_stats(args.reference_profiles_tsv)
            _write_tsv_gz(schema_out, ["feature_name"], schema_rows)
            n_schema_rows = len(schema_rows)
            _write_tsv_gz(stats_out, ["feature_name", "center", "scale"], stats_rows)
        if args.feature_stats_tsv:
            n_stats_rows, _ = _copy_to_tsv_gz(args.feature_stats_tsv, stats_out)
        elif args.feature_schema_tsv:
            _schema_rows, stats_rows = _derive_feature_schema_and_stats(args.reference_profiles_tsv)
            _write_tsv_gz(stats_out, ["feature_name", "center", "scale"], stats_rows)
            n_stats_rows = len(stats_rows)
        else:
            n_stats_rows = n_schema_rows

        bundle_manifest = {
            "bundle_id": bundle_id,
            "assay": "calorimetry",
            "organism": args.organism,
            "summary": {
                "n_reference_profiles": n_profiles,
                "n_reference_metadata_rows": n_metadata,
                "n_feature_schema_rows": n_schema_rows,
                "n_feature_stats_rows": n_stats_rows,
                "n_term_template_rows": n_template_rows,
                "n_phenotype_gene_edges": n_edge_rows,
                "n_term_hierarchy_rows": n_hierarchy_rows,
                "n_mouse_human_ortholog_rows": n_ortholog_rows,
                "term_templates_source": "explicit" if args.term_templates_tsv else "packaged_default",
                "phenotype_gene_edges_source": "explicit" if args.phenotype_gene_edges_tsv else "packaged_default",
                "term_hierarchy_source": "explicit" if args.term_hierarchy_tsv else ("packaged_default" if term_hierarchy_src else "absent"),
                "mouse_human_orthologs_source": "explicit" if getattr(args, "mouse_human_orthologs_tsv", None) else ("packaged_default" if ortholog_src else "absent"),
            },
            "files": {
                "reference_profiles": profiles_out.name,
                "reference_metadata": metadata_out.name,
                "feature_schema": schema_out.name,
                "feature_stats": stats_out.name,
                "term_templates": templates_out.name,
                "phenotype_gene_edges": edges_out.name,
                **({"term_hierarchy": hierarchy_out.name} if term_hierarchy_src else {}),
                **({"mouse_human_orthologs": ortholog_out.name} if ortholog_src else {}),
            },
        }
        bundle_manifest_path = out_dir / f"{bundle_id}.bundle.json"
        write_bundle_manifest(bundle_manifest_path, bundle_manifest)
        bundle_summary = {
            "bundle_id": bundle_id,
            "organism": args.organism,
            "n_reference_profiles": n_profiles,
            "n_reference_metadata_rows": n_metadata,
            "n_feature_schema_rows": n_schema_rows,
            "n_feature_stats_rows": n_stats_rows,
            "n_term_template_rows": n_template_rows,
            "n_phenotype_gene_edges": n_edge_rows,
            "n_term_hierarchy_rows": n_hierarchy_rows,
            "n_mouse_human_ortholog_rows": n_ortholog_rows,
            "term_templates_source": bundle_manifest["summary"]["term_templates_source"],
            "phenotype_gene_edges_source": bundle_manifest["summary"]["phenotype_gene_edges_source"],
            "term_hierarchy_source": bundle_manifest["summary"]["term_hierarchy_source"],
            "mouse_human_orthologs_source": bundle_manifest["summary"]["mouse_human_orthologs_source"],
        }
        (out_dir / "bundle_summary.json").write_text(json.dumps(bundle_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (out_dir / "bundle_summary.txt").write_text(
            "\n".join(
                [
                    f"bundle_id: {bundle_id}",
                    f"organism: {args.organism}",
                    f"n_reference_profiles: {n_profiles}",
                    f"n_reference_metadata_rows: {n_metadata}",
                    f"n_feature_schema_rows: {n_schema_rows}",
                    f"n_feature_stats_rows: {n_stats_rows}",
                    f"n_term_template_rows: {n_template_rows}",
                    f"n_phenotype_gene_edges: {n_edge_rows}",
                    f"n_term_hierarchy_rows: {n_hierarchy_rows}",
                    f"n_mouse_human_ortholog_rows: {n_ortholog_rows}",
                    f"term_templates_source: {bundle_summary['term_templates_source']}",
                    f"phenotype_gene_edges_source: {bundle_summary['phenotype_gene_edges_source']}",
                    f"term_hierarchy_source: {bundle_summary['term_hierarchy_source']}",
                    f"mouse_human_orthologs_source: {bundle_summary['mouse_human_orthologs_source']}",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        distribution_info = None
        if bool(args.write_distribution_artifact):
            included = [
                profiles_out,
                metadata_out,
                schema_out,
                stats_out,
                templates_out,
                edges_out,
                bundle_manifest_path,
                out_dir / "bundle_summary.json",
                out_dir / "bundle_summary.txt",
            ]
            if term_hierarchy_src:
                included.append(hierarchy_out)
            if ortholog_src:
                included.append(ortholog_out)
            distribution_info = _write_distribution_artifact(
                dist_dir=Path(args.distribution_dir).expanduser() if args.distribution_dir else (out_dir / "dist"),
                bundle_id=bundle_id,
                included_files=included,
            )
        return {
            "bundle_id": bundle_id,
            "bundle_manifest": str(bundle_manifest_path),
            "n_profiles": n_profiles,
            "out_dir": str(out_dir),
            **(distribution_info or {}),
        }
