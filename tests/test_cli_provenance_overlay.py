import subprocess, sys, os
from pathlib import Path
DIG = Path(__file__).resolve().parents[1]

def test_cli_provenance_overlay_writes_files(tmp_path):
    env = {**os.environ, "PYTHONPATH": str(DIG / "src")}
    out = tmp_path / "model_out"
    cmd = [sys.executable, "-m", "geneset_extractors.cli", "provenance", "overlay",
           "--pdc_file_manifest_tsv", str(DIG / "tests/data/toy_pdc_file_manifest.tsv"),
           "--prepared_dir", "/prep", "--out_dir", str(out),
           "--operation_script_url", "http://runner"]
    subprocess.run(cmd, cwd=str(DIG), env=env, check=True)
    assert (out / "provenance_overlay.json").exists()
    assert (out / "local_input_source_map.tsv").exists()
