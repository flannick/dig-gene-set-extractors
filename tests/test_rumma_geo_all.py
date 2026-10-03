from pathlib import Path


def test_all_model_provenance_records_consumed_selection_manifest() -> None:
    source = (Path(__file__).resolve().parents[1] / "src/geneset_extractors/extractors/converters/rumma_geo_all.py").read_text(encoding="utf-8")
    assert '"recorded_selection_manifest":str(selection/"selection_manifest.tsv")' in source
    assert '"recorded_selection_manifest":str(query)' not in source
