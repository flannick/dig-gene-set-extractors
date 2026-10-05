"""Pinned IMPC DR18 direct-assertion phenotype gene-set reconstruction."""
from __future__ import annotations

import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context

REQUIRED_COLUMNS = {"marker_symbol", "mp_term_id", "mp_term_name"}
ACRONYMS = {"cd4": "CD4", "cd8": "CD8", "cd25": "CD25", "hdl": "HDL", "ige": "IgE", "igg1": "IgG1", "igg2b": "IgG2b", "klrg1": "KLRG1", "ldl": "LDL", "ly6c": "Ly6C", "nk": "NK", "pq": "PQ", "pr": "PR", "qrs": "QRS", "qt": "QT", "rr": "RR", "st": "ST"}


def _format_term(name: str, mp_id: str) -> str:
    titled = name.strip().title()
    for lower, canonical in ACRONYMS.items():
        titled = re.sub(rf"\b{re.escape(lower)}\b", canonical, titled, flags=re.IGNORECASE)
    return f"{titled} ({mp_id.strip()})"


def _symbol_mapping(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2 and row[0].strip() and row[1].strip():
                mapping[row[0].strip().upper()] = row[1].strip().upper()
    if not mapping:
        raise ValueError(f"no symbol mappings found in {path}")
    return mapping


def _open_csv(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")


def run(args) -> dict[str, object]:
    activate_runtime_context("impc_hz1", getattr(args, "provenance_overlay_json", None))
    assertions, symbol_mapping = Path(args.assertions).resolve(), Path(args.symbol_mapping).resolve()
    if not assertions.is_file(): raise FileNotFoundError(assertions)
    if not symbol_mapping.is_file(): raise FileNotFoundError(symbol_mapping)
    mapping = _symbol_mapping(symbol_mapping); sets: defaultdict[str, set[str]] = defaultdict(set)
    rows_read = rows_accepted = rows_unmapped = 0
    with _open_csv(assertions) as handle:
        reader = csv.DictReader(handle); missing = REQUIRED_COLUMNS.difference(reader.fieldnames or ())
        if missing: raise ValueError(f"assertions missing required columns: {sorted(missing)}")
        for row in reader:
            rows_read += 1
            marker, mp_id, name = (row.get("marker_symbol") or "").strip(), (row.get("mp_term_id") or "").strip(), (row.get("mp_term_name") or "").strip()
            if not marker or not mp_id or not name: continue
            gene = mapping.get(marker.upper())
            if not gene:
                rows_unmapped += 1; continue
            sets[_format_term(name, mp_id)].add(gene); rows_accepted += 1
    emitted = {term: genes for term, genes in sets.items() if len(genes) >= args.min_genes}
    out_dir = Path(args.out_dir).resolve(); out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "genesets.gmt").open("w", encoding="utf-8", newline="\n") as gmt, (out_dir / "geneset.tsv").open("w", encoding="utf-8", newline="") as table:
        writer = csv.DictWriter(table, fieldnames=["term", "gene_symbol", "rank"], delimiter="\t", lineterminator="\n"); writer.writeheader()
        for term in sorted(emitted):
            genes = sorted(emitted[term]); gmt.write("\t".join([term, args.gmt_description, *genes]) + "\n")
            writer.writerows({"term": term, "gene_symbol": gene, "rank": rank} for rank, gene in enumerate(genes, 1))
    summary = {"model_id": args.model_id, "input_release": "IMPC Data Release 18.0", "rows_read": rows_read, "rows_accepted": rows_accepted, "rows_unmapped": rows_unmapped, "terms_before_min_genes": len(sets), "n_gene_sets": len(emitted), "n_memberships": sum(map(len, emitted.values())), "n_genes": len({gene for genes in emitted.values() for gene in genes}), "min_genes": args.min_genes}
    (out_dir / "reconstruction_diagnostics.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    metadata = make_metadata(converter_name="impc_hz1", parameters={"model_id": args.model_id, "source_release": "IMPC Data Release 18.0", "min_genes": args.min_genes, "algorithm": "direct genotype-phenotype assertions; case-normalized Harmonizome mappingFile_2017 symbols; deduplicate term-gene memberships; retain sets with minimum unique genes"}, data_type="phenotype_association", assay="mouse_knockout_phenotype", organism="human", genome_build=args.genome_build, files=[input_file_record(assertions, "impc_dr18_assertions"), input_file_record(symbol_mapping, "harmonizome_mappingfile_2017")], gene_annotation={"mode": "Harmonizome_mappingFile_2017", "gene_id_field": "marker_symbol", "normalization": "case_normalized_lookup_then_uppercase"}, weights={"weight_type": "unweighted", "normalization": {"method": "none"}}, summary=summary, output_files=[{"path": "genesets.gmt", "role": "gmt_library"}, {"path": "geneset.tsv", "role": "selected_program"}, {"path": "reconstruction_diagnostics.json", "role": "reconstruction_diagnostics"}], gmt={"written": True, "path": "genesets.gmt", "prefer_symbol": True, "min_genes": args.min_genes, "max_genes": None, "plans": [{"name": "impc_dr18_hz1_reconstruction", "method": "direct_mp_term_assertion_aggregation", "parameters": {"deterministic_sort": "term then gene"}, "n_genes_emitted": summary["n_memberships"], "token_type": "gene_symbol"}]}, gene_set_description="IMPC Data Release 18 direct mouse knockout phenotype associations normalized with Harmonizome mappingFile_2017.")
    write_metadata(out_dir / "geneset.meta.json", metadata)
    return {**summary, "out_dir": str(out_dir)}
