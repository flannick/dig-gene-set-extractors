"""Finalize streamed LINCS L1000 CP GMT output without loading a term-gene table."""
from __future__ import annotations

import csv
from pathlib import Path

from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context


def _term_direction(set_name: str) -> tuple[str, str, float]:
    if set_name.endswith(" up"):
        return set_name[:-3], "up", 1.0
    if set_name.endswith(" down"):
        return set_name[:-5], "dn", -1.0
    raise ValueError(f"LINCS CP GMT set name must end in ' up' or ' down': {set_name}")


def run(args) -> dict[str, object]:
    activate_runtime_context("lincs_l1000_cp_gmt", getattr(args, "provenance_overlay_json", None))
    source_gmt = Path(args.gmt).resolve()
    if not source_gmt.is_file():
        raise FileNotFoundError(f"Missing streamed LINCS CP GMT: {source_gmt}")
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    description = str(args.gmt_description).strip()
    fields = ["term", "direction", "gene_id", "gene_symbol", "score", "signed_score", "sign", "rank"]
    n_rows = n_sets = 0
    genes: set[str] = set()
    with source_gmt.open(encoding="utf-8", newline="") as source, (out_dir / "genesets.gmt").open("w", encoding="utf-8", newline="\n") as gmt_out, (out_dir / "geneset.tsv").open("w", encoding="utf-8", newline="") as selected_out, (out_dir / "geneset.full.tsv").open("w", encoding="utf-8", newline="") as full_out, (out_dir / "signature_summary.tsv").open("w", encoding="utf-8", newline="") as summary_out:
        selected_writer = csv.DictWriter(selected_out, delimiter="\t", fieldnames=fields, lineterminator="\n")
        full_writer = csv.DictWriter(full_out, delimiter="\t", fieldnames=fields, lineterminator="\n")
        summary_writer = csv.DictWriter(summary_out, delimiter="\t", fieldnames=["source_gmt", "set_name", "description", "gene_count"], lineterminator="\n")
        selected_writer.writeheader()
        full_writer.writeheader()
        summary_writer.writeheader()
        for line_number, line in enumerate(source, start=1):
            values = line.rstrip("\n").split("\t")
            if len(values) < 3:
                raise ValueError(f"Invalid GMT record at {source_gmt}:{line_number}")
            set_name, raw_genes = values[0], [gene for gene in values[2:] if gene]
            term, direction, sign = _term_direction(set_name)
            gmt_out.write("\t".join([set_name, description, *raw_genes]) + "\n")
            summary_writer.writerow({"source_gmt": "genesets.gmt", "set_name": set_name, "description": description, "gene_count": len(raw_genes)})
            n_sets += 1
            for set_rank, gene in enumerate(raw_genes, start=1):
                row = {"term": term, "direction": direction, "gene_id": gene, "gene_symbol": gene, "score": len(raw_genes) - set_rank + 1, "signed_score": sign * (len(raw_genes) - set_rank + 1), "sign": sign, "rank": n_rows + 1}
                selected_writer.writerow(row)
                full_writer.writerow(row)
                n_rows += 1
                genes.add(gene)
    metadata = make_metadata(
        converter_name="lincs_l1000_cp_gmt",
        parameters={"signature_name": args.signature_name, "gmt_description": description, "genes_per_set": args.genes_per_set, "streaming": True},
        data_type="transcriptomics",
        assay="bulk",
        organism=args.organism,
        genome_build=args.genome_build,
        files=[input_file_record(source_gmt, "streamed_lincs_cp_gmt")],
        gene_annotation={"mode": "none", "source": "GCTX row IDs", "gene_id_field": "gene_symbol"},
        weights={"weight_type": "signed_rank", "normalization": {"method": "none", "target_sum": None}, "aggregation": "streamed_per_signature_rank"},
        summary={"n_input_features": n_rows, "n_genes": len(genes), "n_features_assigned": n_rows, "fraction_features_assigned": 1.0, "n_sets_emitted": n_sets},
        output_files=[{"path": "genesets.gmt", "role": "gmt_library"}, {"path": "geneset.tsv", "role": "selected_program"}, {"path": "geneset.full.tsv", "role": "full_scores"}, {"path": "signature_summary.tsv", "role": "signature_summary"}, {"path": "geneset.meta.json", "role": "metadata_json"}],
        gmt={"written": True, "path": "genesets.gmt", "prefer_symbol": True, "min_genes": args.genes_per_set, "max_genes": args.genes_per_set, "plans": [{"name": "streamed_per_signature", "method": "direct_gctx_rank", "parameters": {"description": description, "format": "classic", "signed": True}, "n_genes_emitted": n_rows, "token_type": "gene_symbol", "n_duplicates_dropped": 0}]},
        gene_set_description=args.signature_name,
        upstream_provenance_graph_path=getattr(args, "upstream_provenance_graph_json", None),
        provenance_mirror_local_prefix=getattr(args, "provenance_mirror_local_prefix", None),
        provenance_mirror_remote_prefix=getattr(args, "provenance_mirror_remote_prefix", None),
    )
    write_metadata(out_dir / "geneset.meta.json", metadata)
    return {"n_peaks": n_rows, "n_genes": len(genes), "n_sets": n_sets, "out_dir": str(out_dir)}
