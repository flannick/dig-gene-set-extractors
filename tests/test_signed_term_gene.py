from __future__ import annotations

import json
from pathlib import Path

import yaml

from geneset_extractors.cli import main


def test_signed_term_gene_uses_signature_name_for_metadata_and_dapper_collection(tmp_path: Path) -> None:
    table = tmp_path / "terms.tsv"
    table.write_text(
        "term\tgene_id\tgene_symbol\tscore\tsign\n"
        "Drug_A\tGENE1\tGENE1\t1\t1\n"
        "Drug_A\tGENE2\tGENE2\t1\t1\n",
        encoding="utf-8",
    )
    out_dir = tmp_path / "out"

    assert main([
        "convert", "signed_term_gene", "--table_tsv", str(table), "--out_dir", str(out_dir),
        "--organism", "human", "--genome_build", "hg38", "--signature_name", "LINCS_L1000_Chem_Pert",
        "--gmt_min_genes", "1",
    ]) == 0

    metadata = json.loads((out_dir / "geneset.meta.json").read_text(encoding="utf-8"))
    assert metadata["gene_set"]["name"] == "LINCS_L1000_Chem_Pert"
    payload = yaml.safe_load((out_dir / "geneset.provenance.dapper.yaml").read_text(encoding="utf-8"))
    collection = payload["gene_set_collections"][0]
    rows = payload["gene_sets"]
    assert collection["name"] == "LINCS L1000 Chem Pert"
    assert collection["members"] == [row["id"] for row in rows]
    assert (out_dir / "genesets.dapper-ids.gmt").is_file()
