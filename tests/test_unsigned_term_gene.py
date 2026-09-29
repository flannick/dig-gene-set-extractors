from __future__ import annotations

import json
from pathlib import Path

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
    assert "name: HuBMAP_ASCTB" in dapper
