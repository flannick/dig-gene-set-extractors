import csv
from pathlib import Path

from geneset_extractors.extractors.proteomics.public_prepare import run_public_prepare


def _write(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def test_cdap_log_ratio_columns_join_bare_aliquot_annotations(tmp_path):
    # CDAP phosphosite report: sample columns carry the " Log Ratio" suffix.
    ptm = tmp_path / "phospho.tsv"
    _write(
        ptm,
        "Phosphosite\tA1 Log Ratio\tA2 Log Ratio\tPeptide\tGene\tOrganism\n"
        "NP_1:s10\t1.5\t0.1\tAAAsBBB\tAKT1\thuman\n"
        "NP_2:t20\t1.2\t0.0\tCCCtDDD\tMTOR\thuman\n",
    )
    # sample_annotations keyed by BARE aliquot id (as the PDC biospecimen API provides).
    anno = tmp_path / "anno.tsv"
    _write(
        anno,
        "sample_id_raw\tcondition\tsample_type\taliquot_submitter_id\n"
        "A1\tcase\tPrimary Tumor\tA1\n"
        "A2\tcontrol\tSolid Tissue Normal\tA2\n",
    )
    out = tmp_path / "prepared"
    run_public_prepare(
        input_mode="cdap_files",
        ptm_report_tsv=str(ptm),
        protein_report_tsv=None,
        sample_design_tsv=None,
        sample_annotations_tsv=str(anno),
        pdc_manifest_tsv=None,
        source_dir=None,
        out_dir=str(out),
        organism="human",
        ptm_type="phospho",
        study_id="s1",
        study_label="S1",
    )
    rows = list(csv.DictReader((out / "sample_metadata.tsv").open(encoding="utf-8"), delimiter="\t"))
    conditions = [r["condition"] for r in rows]
    assert "case" in conditions, f"expected a case sample, got {conditions}"
    assert "control" in conditions, f"expected a control sample, got {conditions}"
