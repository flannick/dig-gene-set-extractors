"""Metabolomics Workbench HZ1 gene--metabolite association reconstruction."""
from __future__ import annotations

import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context

REQUIRED_COLUMNS = {"Gene", "Gene ID", "Metabolite", "Metabolite ID"}


def _open_table(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")


def run(args) -> dict[str, object]:
    """Group a Harmonizome MW edge list into deterministic metabolite gene sets."""
    activate_runtime_context("metabolomics_workbench_hz1", getattr(args, "provenance_overlay_json", None))
    edges = Path(args.edges).resolve()
    if not edges.is_file():
        raise FileNotFoundError(edges)
    groups: defaultdict[str, set[str]] = defaultdict(set)
    rows_read = rows_accepted = rows_missing_required_value = 0
    metabolite_ids: dict[str, set[str]] = defaultdict(set)
    with _open_table(edges) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing_columns = REQUIRED_COLUMNS.difference(reader.fieldnames or ())
        if missing_columns:
            raise ValueError(f"edge list missing required columns: {sorted(missing_columns)}")
        for row in reader:
            rows_read += 1
            gene = (row.get("Gene") or "").strip()
            metabolite = (row.get("Metabolite") or "").strip()
            metabolite_id = (row.get("Metabolite ID") or "").strip()
            if not gene or not metabolite or not metabolite_id:
                rows_missing_required_value += 1
                continue
            # This published edge list is already Harmonizome gene-symbol harmonized.
            groups[metabolite].add(gene)
            metabolite_ids[metabolite].add(metabolite_id)
            rows_accepted += 1
    emitted = {term: genes for term, genes in groups.items() if len(genes) >= args.min_genes}
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "genesets.gmt").open("w", encoding="utf-8", newline="\n") as gmt, (out_dir / "geneset.tsv").open("w", encoding="utf-8", newline="") as table:
        writer = csv.DictWriter(table, fieldnames=["metabolite", "metabolite_id", "gene_symbol", "rank"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for metabolite in sorted(emitted):
            genes = sorted(emitted[metabolite])
            gmt.write("\t".join([metabolite, args.gmt_description, *genes]) + "\n")
            writer.writerows(
                {"metabolite": metabolite, "metabolite_id": ";".join(sorted(metabolite_ids[metabolite])), "gene_symbol": gene, "rank": rank}
                for rank, gene in enumerate(genes, 1)
            )
    summary = {
        "model_id": args.model_id, "source_mode": "already_harmonized_harmonizome_edge_list", "rows_read": rows_read,
        "rows_accepted": rows_accepted, "rows_missing_required_value": rows_missing_required_value,
        "terms_before_min_genes": len(groups), "n_gene_sets": len(emitted),
        "n_memberships": sum(len(genes) for genes in emitted.values()),
        "n_genes": len({gene for genes in emitted.values() for gene in genes}), "min_genes": args.min_genes,
    }
    (out_dir / "reconstruction_diagnostics.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    metadata = make_metadata(
        converter_name="metabolomics_workbench_hz1",
        parameters={"model_id": args.model_id, "min_genes": args.min_genes, "algorithm": "deduplicate already-harmonized Harmonizome Gene x Metabolite edges; group by source metabolite label; retain minimum unique-gene sets; sort metabolite and gene axes"},
        data_type="metabolite_enzyme_association", assay="metabolomics_metabolite_enzyme_association", organism="human", genome_build=args.genome_build,
        files=[input_file_record(edges, "harmonizome_mwmetabolites_gene_attribute_edges")],
        gene_annotation={"mode": "source_already_harmonized", "gene_id_field": "Gene ID", "gene_symbol_field": "Gene", "normalization": "none; do not apply a second historical mapping"},
        weights={"weight_type": "unweighted", "normalization": {"method": "none"}}, summary=summary,
        output_files=[{"path": "genesets.gmt", "role": "gmt_library"}, {"path": "geneset.tsv", "role": "selected_program"}, {"path": "reconstruction_diagnostics.json", "role": "reconstruction_diagnostics"}],
        gmt={"written": True, "path": "genesets.gmt", "prefer_symbol": True, "min_genes": args.min_genes, "max_genes": None, "plans": [{"name": "metabolomics_workbench_hz1", "method": "metabolite_gene_edge_grouping", "parameters": {"deterministic_sort": "metabolite then gene"}, "n_genes_emitted": summary["n_memberships"], "token_type": "gene_symbol"}]},
        gene_set_description="Human genes associated with Metabolomics Workbench metabolites, reconstructed from the official Harmonizome processed edge list.",
    )
    write_metadata(out_dir / "geneset.meta.json", metadata)
    return {**summary, "out_dir": str(out_dir)}
