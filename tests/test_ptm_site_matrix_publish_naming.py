import subprocess, sys, os
from pathlib import Path

DIG = Path(__file__).resolve().parents[1]

def _run_convert(out_dir, extra):
    env = {**os.environ, "PYTHONPATH": str(DIG / "src")}
    cmd = [sys.executable, "-m", "geneset_extractors.cli", "convert", "ptm_site_matrix",
           "--ptm_matrix_tsv", str(DIG / "tests/data/toy_ptm_matrix.tsv"),
           "--sample_metadata_tsv", str(DIG / "tests/data/toy_ptm_sample_metadata.tsv"),
           "--protein_matrix_tsv", str(DIG / "tests/data/toy_protein_matrix_gene_level.tsv"),
           "--study_contrast", "condition_a_vs_b", "--condition_a", "case", "--condition_b", "control",
           "--protein_adjustment_run_mode", "compare_if_protein", "--protein_accession_column", "gene_symbol",
           "--select", "top_k", "--top_k", "50", "--gene_aggregation", "signed_topk_mean",
           "--organism", "human", "--genome_build", "human", "--ptm_type", "phospho",
           "--emit_small_gene_sets", "true", "--use_reference_bundle", "false",
           "--out_dir", str(out_dir), *extra]
    subprocess.run(cmd, cwd=str(DIG), env=env, check=True)

def _gmt_names(out_dir):
    names = set()
    for gmt in Path(out_dir).rglob("genesets.gmt"):
        for line in gmt.read_text().splitlines():
            if line.strip():
                names.add(line.split("\t", 1)[0])
    return names

def test_publish_naming_emits_clean_signed_variant_names(tmp_path):
    _run_convert(tmp_path, ["--signature_name", "CPTAC_ClearCellRCC",
                            "--gmt_name_style", "publish", "--gmt_signed_labels", "up,dn"])
    names = _gmt_names(tmp_path)
    expected = {
        "CPTAC_ClearCellRCC_ProteinAdjusted_up", "CPTAC_ClearCellRCC_ProteinAdjusted_dn",
        "CPTAC_ClearCellRCC_Unadjusted_up", "CPTAC_ClearCellRCC_Unadjusted_dn",
    }
    assert expected.issubset(names), f"missing publish names; got {sorted(names)}"
    assert not any("__signature=" in n for n in names), f"verbose scaffold leaked: {sorted(names)}"

def test_verbose_naming_unchanged_by_default(tmp_path):
    _run_convert(tmp_path, ["--signature_name", "ptm_matrix"])
    names = _gmt_names(tmp_path)
    assert any(n.startswith("ptm_site_matrix__signature=") for n in names), sorted(names)

def test_publish_naming_disambiguates_topk_site_comparison_variants(tmp_path):
    _run_convert(tmp_path, ["--signature_name", "CPTAC_Test",
                            "--gmt_name_style", "publish", "--gmt_signed_labels", "up,dn",
                            "--emit_gene_topk_site_comparison", "true",
                            "--gene_topk_site_compare_to", "1"])
    names = _gmt_names(tmp_path)
    expected = {
        "CPTAC_Test_ProteinAdjusted_sites1_up", "CPTAC_Test_ProteinAdjusted_sites1_dn",
        "CPTAC_Test_ProteinAdjusted_sites3_up", "CPTAC_Test_ProteinAdjusted_sites3_dn",
        "CPTAC_Test_Unadjusted_sites1_up", "CPTAC_Test_Unadjusted_sites1_dn",
        "CPTAC_Test_Unadjusted_sites3_up", "CPTAC_Test_Unadjusted_sites3_dn",
    }
    assert names == expected, f"expected 4 distinct sites-disambiguated publish names (x2 signed); got {sorted(names)}"

def test_publish_naming_single_site_cap_has_no_sites_suffix(tmp_path):
    _run_convert(tmp_path, ["--signature_name", "CPTAC_Test",
                            "--gmt_name_style", "publish", "--gmt_signed_labels", "up,dn"])
    names = _gmt_names(tmp_path)
    expected = {
        "CPTAC_Test_ProteinAdjusted_up", "CPTAC_Test_ProteinAdjusted_dn",
        "CPTAC_Test_Unadjusted_up", "CPTAC_Test_Unadjusted_dn",
    }
    assert names == expected, f"single-site-cap publish names must stay unchanged (no _sites token); got {sorted(names)}"
    assert not any("_sites" in n for n in names)
