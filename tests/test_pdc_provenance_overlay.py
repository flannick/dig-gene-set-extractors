from pathlib import Path
from geneset_extractors.extractors.proteomics import pdc_provenance as pdc

DIG = Path(__file__).resolve().parents[1]

def test_build_overlay_maps_roles_to_prepared_paths_with_crdc_ids():
    rows = pdc.read_manifest_tsv(DIG / "tests/data/toy_pdc_file_manifest.tsv")
    ov = pdc.build_overlay(manifest_rows=rows, prepared_dir="/prep", operation_meta={"script_url": "http://x"})
    # Prepared matrices are plain intermediates now: no input decoration, so the
    # workflow-output node and the convert-input node coalesce during the merge.
    assert ov["inputs"] == {}
    # Study-level provenance is still carried on the operation + gene_set nodes.
    assert ov["operation"]["dcc_url"] == "https://pdc.cancer.gov/pdc/study/PDC000128"
    assert ov["gene_set"]["dcc_url"] == "https://pdc.cancer.gov/pdc/study/PDC000128"
    assert ov["operation"]["script_url"] == "http://x"
    assert ov["gene_set"]["drc_url"] == "https://datacommons.cancer.gov/repository/proteomic-data-commons"

def test_write_overlay_emits_overlay_and_source_map(tmp_path):
    rows = pdc.read_manifest_tsv(DIG / "tests/data/toy_pdc_file_manifest.tsv")
    out = pdc.write_overlay(manifest_rows=rows, prepared_dir="/prep", operation_meta={}, out_dir=tmp_path)
    assert out["overlay_json"].exists() and out["source_map_tsv"].exists()
    assert "drs://dg.4DFC/uuid-phos-1" in out["source_map_tsv"].read_text()
