from __future__ import annotations

import json
from pathlib import Path

import yaml

from geneset_extractors.cli import main


def test_unsigned_term_gene_uses_signature_name_for_metadata_and_dapper_collection(tmp_path: Path) -> None:
    table = tmp_path / "terms.tsv"
    table.write_text(
        "term\tgene_id\tgene_symbol\tscore\n"
        "Podocyte\tNPHS1\tNPHS1\t1\n"
        "Podocyte\tWT1\tWT1\t1\n",
        encoding="utf-8",
    )
    out_dir = tmp_path / "out"

    assert main([
        "convert", "unsigned_term_gene", "--table_tsv", str(table), "--out_dir", str(out_dir),
        "--organism", "human", "--genome_build", "hg38", "--signature_name", "HuBMAP_ASCTB",
        "--gmt_min_genes", "1",
    ]) == 0

    metadata = json.loads((out_dir / "geneset.meta.json").read_text(encoding="utf-8"))
    assert metadata["gene_set"]["name"] == "HuBMAP_ASCTB"
    dapper = (out_dir / "geneset.provenance.dapper.yaml").read_text(encoding="utf-8")
    payload = yaml.safe_load(dapper)
    collection = payload["gene_set_collections"][0]
    rows = payload["gene_sets"]
    companion = out_dir / "genesets.dapper-ids.gmt"
    assert collection["name"] == "HuBMAP ASCTB"
    assert collection["members"] == [row["id"] for row in rows]
    assert companion.is_file()
    companion_file = next(node for node in payload["files"] if node["filename"] == companion.name)
    assert collection["has_gmt_file"] == companion_file["id"]
    for row, line in zip(rows, companion.read_text(encoding="utf-8").splitlines(), strict=True):
        assert row["id"] == row["gmt_entry"] == line.split("\t", 1)[0]
        assert row["in_gmt_file"] == companion_file["id"]
        assert row["in_gene_set_collection"] == [collection["id"]]
