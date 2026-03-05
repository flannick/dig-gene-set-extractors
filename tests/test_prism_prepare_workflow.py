from pathlib import Path

from geneset_extractors.workflows.prism_prepare import run as run_prism_prepare


class Args:
    out_dir = "tests/tmp/prism_prepare"
    release = "test"
    matrix_url = "tests/data/toy_prism_matrix.csv"
    treatment_info_url = "tests/data/toy_prism_treatment_info.csv"
    cell_line_info_url = "tests/data/toy_prism_cell_line_info.csv"
    matrix_file_id = ""
    treatment_info_file_id = ""
    cell_line_info_file_id = ""
    depmap_file_id_url_template = "https://depmap.org/portal/download/api/downloads?file_name={file_id}"
    max_retries = 2
    retry_backoff_sec = 0.01
    user_agent = "geneset-extractors-test"


def test_prism_prepare_local_sources(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "prism_prepare")
    result = run_prism_prepare(args)
    assert int(result["n_response_rows"]) > 0
    out_dir = Path(result["out_dir"])
    assert (out_dir / "response_long.tsv").exists()
    assert (out_dir / "groups.tsv").exists()
    assert (out_dir / "drug_targets.tsv").exists()
    assert (out_dir / "prepare_summary.json").exists()

