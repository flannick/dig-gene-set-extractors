"""Deterministic GlyGen gene-set reconstructions from pinned source snapshots."""
from __future__ import annotations

import csv
import hashlib
import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context


API_TEMPLATE = "https://api.glygen.org/glycan/detail/{accession}/"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return "sha256:" + digest.hexdigest()


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_outputs(out_dir: Path, gene_sets: dict[str, set[str]], description: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "genesets.gmt").open("w", encoding="utf-8", newline="\n") as gmt, (out_dir / "geneset.tsv").open("w", encoding="utf-8", newline="") as table:
        writer = csv.DictWriter(table, fieldnames=["term", "gene_symbol", "rank"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for term in sorted(gene_sets):
            genes = sorted(gene_sets[term])
            gmt.write("\t".join([term, description, *genes]) + "\n")
            for rank, gene in enumerate(genes, start=1):
                writer.writerow({"term": term, "gene_symbol": gene, "rank": rank})


def _summary(gene_sets: dict[str, set[str]]) -> dict[str, int]:
    return {"n_gene_sets": len(gene_sets), "n_memberships": sum(map(len, gene_sets.values())), "n_genes": len(set().union(*gene_sets.values())) if gene_sets else 0}


def run_glycosylated_proteins(args) -> dict[str, object]:
    """Reconstruct glycan-associated proteins from the v1.12.1 CSV release."""
    activate_runtime_context("glygen_glycosylated_proteins", getattr(args, "provenance_overlay_json", None))
    citation_paths = [Path(getattr(args, name)).resolve() for name in ("unicarbkb", "harvard", "glyconnect")]
    masterlist = Path(args.masterlist).resolve()
    for path in [*citation_paths, masterlist]:
        if not path.is_file():
            raise FileNotFoundError(path)
    pairs: set[tuple[str, str]] = set()
    for path in citation_paths:
        for row in _read_rows(path):
            protein, glycan = row.get("uniprotkb_canonical_ac", "").strip(), row.get("glytoucan_ac", "").strip()
            if protein and glycan:
                pairs.add((protein, glycan))
    mapping: dict[str, str] = {}
    for row in _read_rows(masterlist):
        protein, gene = row.get("uniprotkb_canonical_ac", "").strip(), row.get("gene_name", "").strip().upper()
        if protein and gene and protein not in mapping:
            mapping[protein] = gene
    grouped: dict[str, set[str]] = defaultdict(set)
    for protein, glycan in pairs:
        if protein in mapping:
            grouped[glycan].add(mapping[protein])
    emitted = {term: genes for term, genes in grouped.items() if len(genes) >= args.min_genes}
    out_dir = Path(args.out_dir).resolve()
    _write_outputs(out_dir, emitted, args.gmt_description)
    diagnostics = {"model_id": "glycosylated_proteins", "release": "GlyGen v1.12.1", "min_genes": args.min_genes, "n_unique_protein_glycan_pairs": len(pairs), "n_sets_before_min_genes": len(grouped), **_summary(emitted)}
    (out_dir / "reconstruction_diagnostics.json").write_text(json.dumps(diagnostics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    files = [input_file_record(path, role) for path, role in zip(citation_paths, ("unicarbkb_glycosylation_citations", "harvard_glycosylation_citations", "glyconnect_glycosylation_citations"))] + [input_file_record(masterlist, "human_protein_masterlist")]
    metadata = make_metadata(converter_name="glygen_glycosylated_proteins", parameters={"model_id": "glycosylated_proteins", "source_release": "v1.12.1", "min_genes": args.min_genes, "algorithm": "deduplicate protein-glycan pairs; map canonical UniProt accessions using same-release masterlist; uppercase symbols; retain sets with minimum unique genes"}, data_type="proteomics", assay="glycosylation_annotation", organism="human", genome_build=args.genome_build, files=files, gene_annotation={"mode": "same_release_GlyGen_masterlist", "gene_id_field": "gene_symbol", "normalization": "strip_uppercase"}, weights={"weight_type": "unweighted", "normalization": {"method": "none"}}, summary=diagnostics, output_files=[{"path": "genesets.gmt", "role": "gmt_library"}, {"path": "geneset.tsv", "role": "selected_program"}, {"path": "reconstruction_diagnostics.json", "role": "reconstruction_diagnostics"}], gmt={"written": True, "path": "genesets.gmt", "prefer_symbol": True, "min_genes": args.min_genes, "max_genes": None, "plans": [{"name": "glygen_v1_12_1_scientific_reimplementation", "method": "glycan_to_glycosylated_proteins", "parameters": {"deterministic_sort": "term then gene"}, "n_genes_emitted": diagnostics["n_memberships"], "token_type": "gene_symbol"}]}, gene_set_description="GlyGen v1.12.1 glycosylated proteins grouped by GlyTouCan accession.")
    write_metadata(out_dir / "geneset.meta.json", metadata)
    return {**diagnostics, "out_dir": str(out_dir)}


def run_glycan_synthesizing_enzymes(args) -> dict[str, object]:
    """Transform cached GlyGen glycan-detail responses; never contacts the API."""
    activate_runtime_context("glygen_glycan_synthesizing_enzymes", getattr(args, "provenance_overlay_json", None))
    manifest = Path(args.cache_manifest).resolve()
    if not manifest.is_file():
        raise FileNotFoundError(manifest)
    gene_sets: dict[str, set[str]] = {}
    cached_files: list[Path] = []
    with manifest.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("status") not in {"cached", "fetched"} or not row.get("cache_file", "").strip():
                continue
            path = Path(row["cache_file"])
            if not path.is_absolute(): path = manifest.parent / path
            if not path.is_file(): raise FileNotFoundError(f"cache manifest references missing response: {path}")
            payload = json.loads(path.read_text(encoding="utf-8"))
            genes = {str(record.get("gene", "")).strip().upper() for record in (payload.get("enzyme") or []) if str(record.get("tax_id", "")) == "9606" and str(record.get("gene", "")).strip()}
            if genes: gene_sets[f"glytoucan:{row['accession'].strip()}"] = genes
            cached_files.append(path)
    out_dir = Path(args.out_dir).resolve()
    _write_outputs(out_dir, gene_sets, args.gmt_description)
    diagnostics = {"model_id": "glycan_synthesizing_enzymes", "source": "GlyGen glycan-detail API / Sandbox biosynthetic-enzyme annotations", "human_tax_id": "9606", "xref_key": "glycan_xref_sandbox", **_summary(gene_sets)}
    (out_dir / "reconstruction_diagnostics.json").write_text(json.dumps(diagnostics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    metadata = make_metadata(converter_name="glygen_glycan_synthesizing_enzymes", parameters={"model_id": "glycan_synthesizing_enzymes", "source": "GlyGen glycan-detail API enzyme[]; xref_key=glycan_xref_sandbox", "tax_id": "9606", "algorithm": "cached response enzyme[] records filtered to human, normalized uppercase, deduplicated within GlyTouCan accession"}, data_type="proteomics", assay="glycan_biosynthesis_annotation", organism="human", genome_build=args.genome_build, files=[input_file_record(manifest, "glygen_api_cache_manifest")] + [input_file_record(path, "glygen_glycan_detail_response") for path in cached_files], gene_annotation={"mode": "GlyGen_API_gene_symbol", "gene_id_field": "gene", "normalization": "strip_uppercase"}, weights={"weight_type": "unweighted", "normalization": {"method": "none"}}, summary=diagnostics, output_files=[{"path": "genesets.gmt", "role": "gmt_library"}, {"path": "geneset.tsv", "role": "selected_program"}, {"path": "reconstruction_diagnostics.json", "role": "reconstruction_diagnostics"}], gmt={"written": True, "path": "genesets.gmt", "prefer_symbol": True, "min_genes": 1, "max_genes": None, "plans": [{"name": "glygen_sandbox_biosynthetic_enzymes", "method": "glycan_to_human_biosynthetic_enzyme", "parameters": {"deterministic_sort": "term then gene"}, "n_genes_emitted": diagnostics["n_memberships"], "token_type": "gene_symbol"}]}, gene_set_description="GlyGen Sandbox biosynthetic-enzyme annotations grouped by GlyTouCan accession.")
    write_metadata(out_dir / "geneset.meta.json", metadata)
    return {**diagnostics, "out_dir": str(out_dir)}


def run_glycan_synthesizing_enzymes_acquire(args) -> dict[str, object]:
    """Acquire/cache API responses from an independently supplied accession snapshot."""
    accessions = Path(args.accessions_tsv).resolve()
    cache_dir, manifest = Path(args.cache_dir).resolve(), Path(args.manifest).resolve()
    if not accessions.is_file(): raise FileNotFoundError(accessions)
    cache_dir.mkdir(parents=True, exist_ok=True); manifest.parent.mkdir(parents=True, exist_ok=True)
    ids = sorted({row[args.accession_column].strip().removeprefix("glytoucan:") for row in _read_rows(accessions) if row.get(args.accession_column, "").strip()})
    rows = []
    for accession in ids:
        cache = cache_dir / f"{accession}.json"; url = API_TEMPLATE.format(accession=urllib.parse.quote(accession))
        try:
            if cache.is_file(): status = "cached"
            else:
                error = None
                for attempt in range(1, args.retries + 1):
                    try:
                        request = urllib.request.Request(url, headers={"User-Agent": "geneset-extractors-glygen/1.0", "Accept": "application/json"})
                        with urllib.request.urlopen(request, timeout=args.timeout_seconds) as response: cache.write_bytes(response.read())
                        error = None; break
                    except Exception as exc:  # preserve the failure in the manifest after all retries
                        error = exc
                        if attempt < args.retries: time.sleep(max(args.pause_seconds, 1.0) * attempt)
                if error is not None: raise error
                status = "fetched"; time.sleep(args.pause_seconds)
            rows.append({"accession": accession, "request_url": url, "acquired_at": datetime.now(timezone.utc).isoformat(), "cache_file": str(cache.relative_to(manifest.parent)) if cache.is_relative_to(manifest.parent) else str(cache), "sha256": _sha256(cache), "status": status, "error": ""})
        except Exception as exc:
            rows.append({"accession": accession, "request_url": url, "acquired_at": datetime.now(timezone.utc).isoformat(), "cache_file": str(cache), "sha256": "", "status": "failed", "error": str(exc)})
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["accession", "request_url", "acquired_at", "cache_file", "sha256", "status", "error"], delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    return {"n_accessions": len(ids), "n_failed": sum(row["status"] == "failed" for row in rows), "manifest": str(manifest)}
