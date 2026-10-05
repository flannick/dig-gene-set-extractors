"""Deterministic IDG Enrichr-library acquisition and conversion."""
from __future__ import annotations

import csv
import hashlib
import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context


ENRICHR_ENDPOINT = "https://maayanlab.cloud/Enrichr/geneSetLibrary"
MODEL_LIBRARIES = {
    "idg_drug_targets_2022": "IDG_Drug_Targets_2022",
    "idg_archs4_coexp": "ARCHS4_IDG_Coexp",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return "sha256:" + digest.hexdigest()


def parse_enrichr_gmt(path: Path) -> dict[str, list[str]]:
    """Parse ordinary GMT, retaining source member order and omitting empty sets."""
    gene_sets: dict[str, list[str]] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        fields = raw_line.split("\t")
        term = fields[0].strip() if fields else ""
        members = [member.strip() for member in fields[2:] if member.strip()]
        if not term or not members:
            continue
        if term in gene_sets:
            raise ValueError(f"duplicate GMT term {term!r} in {path}")
        gene_sets[term] = members
    return gene_sets


def _write_outputs(out_dir: Path, gene_sets: dict[str, list[str]], description: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "genesets.gmt").open("w", encoding="utf-8", newline="\n") as gmt, (out_dir / "geneset.tsv").open("w", encoding="utf-8", newline="") as table:
        writer = csv.DictWriter(table, fieldnames=["term", "gene_symbol", "rank"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for term in sorted(gene_sets):
            genes = gene_sets[term]
            gmt.write("\t".join([term, description, *genes]) + "\n")
            writer.writerows({"term": term, "gene_symbol": gene, "rank": rank} for rank, gene in enumerate(genes, start=1))


def _acquire(model: str, destination: Path, timeout_seconds: int) -> tuple[Path, str]:
    library = MODEL_LIBRARIES[model]
    url = f"{ENRICHR_ENDPOINT}?{urllib.parse.urlencode({'mode': 'text', 'libraryName': library})}"
    request = urllib.request.Request(url, headers={"User-Agent": "geneset-extractors-idg/1.0", "Accept": "text/plain"})
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
        payload = response.read()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(payload)
    acquisition = {
        "model": model,
        "library_name": library,
        "requested_url": url,
        "acquired_at": datetime.now(timezone.utc).isoformat(),
        "artifact": str(destination),
        "sha256": _sha256(destination),
    }
    destination.with_suffix(destination.suffix + ".acquisition.json").write_text(json.dumps(acquisition, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return destination, url


def run(args) -> dict[str, object]:
    model = getattr(args, "model", None) or args.converter
    if model not in MODEL_LIBRARIES:
        raise ValueError(f"unsupported IDG model: {model}")
    activate_runtime_context(model, getattr(args, "provenance_overlay_json", None))
    out_dir = Path(args.out_dir).resolve()
    input_value = getattr(args, "input_gmt", None)
    if input_value:
        input_gmt = Path(input_value).resolve()
        if not input_gmt.is_file():
            raise FileNotFoundError(input_gmt)
        source_url = getattr(args, "source_url", None) or input_gmt.as_uri()
    else:
        input_gmt, source_url = _acquire(model, out_dir / "acquisition" / f"{MODEL_LIBRARIES[model]}.gmt", args.timeout_seconds)
    gene_sets = parse_enrichr_gmt(input_gmt)
    _write_outputs(out_dir, gene_sets, args.gmt_description)
    summary = {
        "model": model,
        "library_name": MODEL_LIBRARIES[model],
        "n_gene_sets": len(gene_sets),
        "n_memberships": sum(len(genes) for genes in gene_sets.values()),
        "n_genes": len({gene for genes in gene_sets.values() for gene in genes}),
        "empty_records_omitted": sum(1 for line in input_gmt.read_text(encoding="utf-8").splitlines() if len(line.split("\t")) >= 2 and not any(member.strip() for member in line.split("\t")[2:])),
    }
    (out_dir / "reconstruction_diagnostics.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    metadata = make_metadata(
        converter_name=model,
        parameters={"model": model, "enrichr_library": MODEL_LIBRARIES[model], "source_url": source_url, "algorithm": "parse named Enrichr GMT; discard records without members; preserve member spelling and order"},
        data_type="knowledgebase",
        assay="drug_target_annotation" if model == "idg_drug_targets_2022" else "gene_coexpression",
        organism="human",
        genome_build=args.genome_build,
        files=[input_file_record(input_gmt, "enrichr_gmt", resource_record={"url": source_url, "version": _sha256(input_gmt)})],
        gene_annotation={"mode": "source_gene_symbols", "gene_id_field": "gene_symbol", "normalization": "none"},
        weights={"weight_type": "unweighted", "normalization": {"method": "none"}},
        summary=summary,
        output_files=[{"path": "genesets.gmt", "role": "gmt_library"}, {"path": "geneset.tsv", "role": "selected_program"}, {"path": "reconstruction_diagnostics.json", "role": "reconstruction_diagnostics"}],
        gmt={"written": True, "path": "genesets.gmt", "prefer_symbol": True, "min_genes": 1, "max_genes": None, "plans": [{"name": "named_enrichr_library", "method": "deterministic_gmt_parse", "parameters": {"omit_empty_sets": True, "term_order": "lexicographic", "member_order": "source"}, "n_genes_emitted": summary["n_memberships"], "token_type": "gene_symbol"}]},
        gene_set_description=f"Enrichr {MODEL_LIBRARIES[model]} gene sets.",
    )
    write_metadata(out_dir / "geneset.meta.json", metadata)
    return {**summary, "out_dir": str(out_dir)}
