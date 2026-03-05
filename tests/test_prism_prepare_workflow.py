import hashlib
import json
from pathlib import Path

from geneset_extractors.workflows.prism_prepare import run as run_prism_prepare, sniff_delimited_payload


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
    subset_seed = 0
    max_cell_lines_total = None
    max_compounds_total = None
    max_cell_lines_per_group = None
    balance_by = None
    min_per_balance_bin = 5
    keep_raw_downloads = False


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


def test_prism_prepare_content_sniffing():
    assert sniff_delimited_payload(b'{"a":1}')["kind"] == "json_like"
    assert sniff_delimited_payload(b"<html><body>x</body></html>")["kind"] == "html_like"
    assert sniff_delimited_payload(b"row_name,COL_A\nCL1,-2.0\n")["kind"] == "delimited_text"


def test_prism_prepare_deterministic_balanced_subset(tmp_path: Path):
    args_a = Args()
    args_a.out_dir = str(tmp_path / "prep_a")
    args_a.subset_seed = 13
    args_a.balance_by = "primary_tissue"
    args_a.min_per_balance_bin = 1
    args_a.max_cell_lines_per_group = 1
    args_a.max_cell_lines_total = 2
    args_a.max_compounds_total = 2
    result_a = run_prism_prepare(args_a)
    out_a = Path(result_a["out_dir"])
    digest_a = hashlib.sha256((out_a / "response_long.tsv").read_bytes()).hexdigest()

    args_b = Args()
    args_b.out_dir = str(tmp_path / "prep_b")
    args_b.subset_seed = 13
    args_b.balance_by = "primary_tissue"
    args_b.min_per_balance_bin = 1
    args_b.max_cell_lines_per_group = 1
    args_b.max_cell_lines_total = 2
    args_b.max_compounds_total = 2
    result_b = run_prism_prepare(args_b)
    out_b = Path(result_b["out_dir"])
    digest_b = hashlib.sha256((out_b / "response_long.tsv").read_bytes()).hexdigest()

    assert digest_a == digest_b
    summary = json.loads((out_a / "prepare_summary.json").read_text(encoding="utf-8"))
    assert int(summary["subset"]["n_samples_selected"]) == 2
