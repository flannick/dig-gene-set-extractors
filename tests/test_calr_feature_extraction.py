import csv
from pathlib import Path

from geneset_extractors.extractors.calorimetry.contrasts import build_group_contrasts
from geneset_extractors.extractors.calorimetry.features import CalorimetryFeatureConfig, SubjectSummary, extract_subject_features
from geneset_extractors.extractors.calorimetry.io import read_calr_data_csv, read_exclusions_tsv, read_session_csv


def test_calr_feature_extraction_selected_window_and_exclusions():
    fieldnames, rows = read_calr_data_csv("tests/data/toy_calr_data.csv")
    session_fieldnames, session_rows = read_session_csv("tests/data/toy_calr_session.csv")
    exclusions = read_exclusions_tsv("tests/data/toy_calr_exclusions.tsv")
    result = extract_subject_features(
        data_rows=rows,
        data_fieldnames=fieldnames,
        session_rows=session_rows,
        session_fieldnames=session_fieldnames,
        exclusion_rows=exclusions,
        cfg=CalorimetryFeatureConfig(
            analysis_start_hour=None,
            analysis_end_hour=None,
            photoperiod_lights_on_hour=None,
            photoperiod_hours_light=12.0,
            exclusions_tsv="tests/data/toy_calr_exclusions.tsv",
            exploratory_without_session=True,
        ),
    )
    assert result.metadata_summary["session_mode"] == "explicit"
    assert result.metadata_summary["analysis_window_source"] == "session"
    assert result.metadata_summary["ee_source"] == "provided"
    assert result.metadata_summary["rer_source"] == "provided"
    assert result.metadata_summary["n_subjects"] == 4
    assert result.metadata_summary["n_rows_selected"] == 14
    assert any("manual_tail_exclusion" in warning for warning in result.warnings)
    assert "vo2_mean_selected" in result.feature_names
    assert "xytot_amp1" in result.feature_names
    assert result.subjects["KO2"].feature_values["feed_auc_selected"] < result.subjects["KO1"].feature_values["feed_auc_selected"]


def test_calr_feature_extraction_exploratory_without_session_warns(capsys):
    fieldnames, rows = read_calr_data_csv("tests/data/toy_calr_data.csv")
    result = extract_subject_features(
        data_rows=rows,
        data_fieldnames=fieldnames,
        session_rows=None,
        session_fieldnames=None,
        exclusion_rows=[],
        cfg=CalorimetryFeatureConfig(
            analysis_start_hour=None,
            analysis_end_hour=None,
            photoperiod_lights_on_hour=6.0,
            photoperiod_hours_light=12.0,
            exclusions_tsv=None,
            exploratory_without_session=True,
        ),
    )
    captured = capsys.readouterr()
    assert result.metadata_summary["session_mode"] == "exploratory"
    assert result.metadata_summary["analysis_window_source"] == "full_trace"
    assert "not fully CalR-equivalent" in captured.err


def test_calr_contrasts_warn_when_mass_covariate_missing(tmp_path: Path, capsys):
    src = Path("tests/data/toy_calr_data.csv")
    dst = tmp_path / "no_mass.csv"
    with src.open("r", encoding="utf-8", newline="") as fh_in, dst.open("w", encoding="utf-8", newline="") as fh_out:
        reader = csv.DictReader(fh_in)
        fieldnames = [name for name in reader.fieldnames if name != "subject.mass"]
        writer = csv.DictWriter(fh_out, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            row.pop("subject.mass", None)
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    session_path = tmp_path / "session_no_mass.csv"
    session_path.write_text(
        "id,group1,analysis_start_hour,analysis_end_hour,lights_on_hour\nWT1,WT,0,18,6\nWT2,WT,0,18,6\nKO1,KO,0,18,6\nKO2,KO,0,18,6\n",
        encoding="utf-8",
    )
    fieldnames, rows = read_calr_data_csv(dst)
    session_fieldnames, session_rows = read_session_csv(session_path)
    extracted = extract_subject_features(
        data_rows=rows,
        data_fieldnames=fieldnames,
        session_rows=session_rows,
        session_fieldnames=session_fieldnames,
        exclusion_rows=[],
        cfg=CalorimetryFeatureConfig(
            analysis_start_hour=None,
            analysis_end_hour=None,
            photoperiod_lights_on_hour=None,
            photoperiod_hours_light=12.0,
            exclusions_tsv=None,
            exploratory_without_session=True,
        ),
    )
    contrasts, warnings = build_group_contrasts(subjects_by_id=extracted.subjects, explicit_mass_covariate=None, min_group_size=2)
    captured = capsys.readouterr()
    assert contrasts
    assert any("no suitable mass/body-composition covariate" in warning for warning in warnings)
    assert "exploratory fallback mode" in captured.err


def test_calr_interaction_kept_vs_dropped():
    def _subjects(values_a, values_b):
        subjects = {}
        for idx, (mass, value) in enumerate(values_a, start=1):
            subjects[f"A{idx}"] = SubjectSummary(
                subject_id=f"A{idx}",
                group="A",
                run_id="RUN1",
                feature_values={"vo2_mean_selected": float(value)},
                metadata={"lean_mass": float(mass)},
            )
        for idx, (mass, value) in enumerate(values_b, start=1):
            subjects[f"B{idx}"] = SubjectSummary(
                subject_id=f"B{idx}",
                group="B",
                run_id="RUN1",
                feature_values={"vo2_mean_selected": float(value)},
                metadata={"lean_mass": float(mass)},
            )
        return subjects

    dropped_subjects = _subjects([(20, 40), (22, 44), (24, 48)], [(20, 42), (22, 46), (24, 50)])
    contrasts_dropped, _warnings = build_group_contrasts(subjects_by_id=dropped_subjects, explicit_mass_covariate="lean_mass", min_group_size=2)
    dropped_model = contrasts_dropped[0].feature_models["vo2_mean_selected"]
    assert dropped_model["interaction_tested"] is True
    assert dropped_model["interaction_retained"] is False

    kept_subjects = _subjects([(20, 40), (22, 44), (24, 48)], [(20, 35), (22, 50), (24, 65)])
    contrasts_kept, _warnings = build_group_contrasts(subjects_by_id=kept_subjects, explicit_mass_covariate="lean_mass", min_group_size=2)
    kept_model = contrasts_kept[0].feature_models["vo2_mean_selected"]
    assert kept_model["interaction_tested"] is True
    assert kept_model["interaction_retained"] is True
