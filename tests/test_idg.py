from argparse import Namespace

from geneset_extractors.extractors.converters.idg import MODEL_LIBRARIES, parse_enrichr_gmt, run


def test_parse_enrichr_gmt_omits_empty_records_and_preserves_members(tmp_path):
    source = tmp_path / "source.gmt"
    source.write_text("DrugA\tdescription\tGENE1\tGENE2\nEmpty\tdescription\n\t\t\nDrugB\tdescription\tGENE3\n")
    assert parse_enrichr_gmt(source) == {"DrugA": ["GENE1", "GENE2"], "DrugB": ["GENE3"]}


def test_idg_conversion_is_deterministic_and_has_correct_identity(tmp_path):
    source = tmp_path / "source.gmt"
    source.write_text("Z\tdescription\tGENE2\tGENE1\nEmpty\tdescription\nA\tdescription\tGENE3\n")
    first, second = tmp_path / "first", tmp_path / "second"
    base = dict(model="idg_archs4_coexp", input_gmt=source, source_url="https://example.test/ARCHS4_IDG_Coexp", genome_build="hg38", gmt_description="test", provenance_overlay_json=None, timeout_seconds=60)
    result = run(Namespace(**base, out_dir=first))
    run(Namespace(**base, out_dir=second))
    assert result["library_name"] == "ARCHS4_IDG_Coexp"
    assert result["empty_records_omitted"] == 1
    assert (first / "genesets.gmt").read_bytes() == (second / "genesets.gmt").read_bytes()
    assert (first / "genesets.gmt").read_text() == "A\ttest\tGENE3\nZ\ttest\tGENE2\tGENE1\n"
    assert MODEL_LIBRARIES == {"idg_drug_targets_2022": "IDG_Drug_Targets_2022", "idg_archs4_coexp": "ARCHS4_IDG_Coexp"}
