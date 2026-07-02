from pathlib import Path
from geneset_extractors.extractors.proteomics import pdc_provenance as pdc

DIG = Path(__file__).resolve().parents[1]

def test_build_overlay_maps_roles_to_prepared_paths_with_crdc_ids():
    rows = pdc.read_manifest_tsv(DIG / "tests/data/toy_pdc_file_manifest.tsv")
    ov = pdc.build_overlay(manifest_rows=rows, prepared_dir="/prep", operation_meta={"script_url": "http://x"})
    ptm = ov["inputs"]["/prep/ptm_matrix.tsv"]
    assert ptm["local_id"] == "drs://dg.4DFC/uuid-phos-1"
    assert ptm["dcc_url"] == "https://pdc.cancer.gov/pdc/study/PDC000128"
    assert ptm["drc_url"] == "https://datacommons.cancer.gov/repository/proteomic-data-commons"
    assert ov["inputs"]["/prep/protein_matrix.tsv"]["persistent_id"] == "uuid-prot-1"
    assert "/prep/sample_metadata.tsv" in ov["inputs"]
    assert ov["operation"]["script_url"] == "http://x"
    assert ov["gene_set"]["dcc_url"].endswith("PDC000128")

def test_write_overlay_emits_overlay_and_source_map(tmp_path):
    rows = pdc.read_manifest_tsv(DIG / "tests/data/toy_pdc_file_manifest.tsv")
    out = pdc.write_overlay(manifest_rows=rows, prepared_dir="/prep", operation_meta={}, out_dir=tmp_path)
    assert out["overlay_json"].exists() and out["source_map_tsv"].exists()
    assert "drs://dg.4DFC/uuid-phos-1" in out["source_map_tsv"].read_text()
