import csv
import json
from pathlib import Path

from geneset_extractors.converters import morphology_profile_query
from tests.test_morphology_profile_query_converter import Args


def _genes(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as fh:
        return [row["gene_id"] for row in csv.DictReader(fh, delimiter="\t")]


def _gene_scores(path: Path) -> dict[str, float]:
    with path.open("r", encoding="utf-8") as fh:
        return {row["gene_id"]: float(row["score"]) for row in csv.DictReader(fh, delimiter="\t")}


def test_morphology_specificity_benchmark_target_recovery_beats_generic_hub(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    refs.write_text(
        "perturbation_id\tf1\tf2\n"
        "HUB1\t1.0\t0.01\n"
        "HUB2\t0.99\t0.02\n"
        "SPEC\t0.80\t0.20\n",
        encoding="utf-8",
    )
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\thub_score\tqc_weight\tis_control\n"
        "HUB1\tcompound\tHUB1\t10.0\t1.0\tfalse\n"
        "HUB2\tcompound\tHUB2\t9.0\t1.0\tfalse\n"
        "SPEC\tcompound\tSPEC\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "compound_id\tgene_symbol\tweight\n"
        "HUB1\tP2RY12\t1.0\n"
        "HUB2\tOPRM1\t1.0\n"
        "SPEC\tKCNN4\t1.0\n",
        encoding="utf-8",
    )
    args = Args()
    args.query_profiles_tsv = str(query)
    args.query_metadata_tsv = None
    args.group_query_by = None
    args.reference_profiles_tsv = str(refs)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "out")
    args.polarity = "similar"
    args.top_k = 3
    args.gmt_topk_list = "3"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    genes = _genes(Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv")
    assert genes[0] == "KCNN4"


def test_morphology_specificity_benchmark_default_suppresses_recurrent_generic_genes(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    ref_lines = ["perturbation_id\tf1\tf2"]
    meta_lines = ["perturbation_id\tperturbation_type\tcompound_id\thub_score\tqc_weight\tis_control"]
    target_lines = ["compound_id\tgene_symbol\tweight"]
    generic_genes = ["P2RY12", "ADA", "OPRM1", "TBXAS1", "ANXA1", "CDC25A"]
    for idx, gene in enumerate(generic_genes, start=1):
        ref_lines.append(f"H{idx}\t1.0\t0.0")
        meta_lines.append(f"H{idx}\tcompound\tH{idx}\t8.0\t1.0\tfalse")
        target_lines.append(f"H{idx}\t{gene}\t1.0")
    ref_lines.append("SPEC\t0.85\t0.15")
    meta_lines.append("SPEC\tcompound\tSPEC\t0.1\t1.0\tfalse")
    target_lines.append("SPEC\tNTRK1\t1.0")
    refs = tmp_path / "refs.tsv"
    refs.write_text("\n".join(ref_lines) + "\n", encoding="utf-8")
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text("\n".join(meta_lines) + "\n", encoding="utf-8")
    targets = tmp_path / "targets.tsv"
    targets.write_text("\n".join(target_lines) + "\n", encoding="utf-8")

    args_default = Args()
    args_default.query_profiles_tsv = str(query)
    args_default.query_metadata_tsv = None
    args_default.group_query_by = None
    args_default.reference_profiles_tsv = str(refs)
    args_default.reference_metadata_tsv = str(ref_meta)
    args_default.compound_targets_tsv = str(targets)
    args_default.feature_stats_tsv = None
    args_default.feature_schema_tsv = None
    args_default.out_dir = str(tmp_path / "default")
    args_default.polarity = "similar"
    args_default.top_k = 7
    args_default.gmt_topk_list = "7"
    args_default.gmt_min_genes = 1
    args_default.emit_small_gene_sets = True
    morphology_profile_query.run(args_default)
    genes_default = _genes(Path(args_default.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv")

    args_broad = Args()
    args_broad.query_profiles_tsv = str(query)
    args_broad.query_metadata_tsv = None
    args_broad.group_query_by = None
    args_broad.reference_profiles_tsv = str(refs)
    args_broad.reference_metadata_tsv = str(ref_meta)
    args_broad.compound_targets_tsv = str(targets)
    args_broad.feature_stats_tsv = None
    args_broad.feature_schema_tsv = None
    args_broad.out_dir = str(tmp_path / "broad")
    args_broad.polarity = "similar"
    args_broad.max_reference_neighbors = 0
    args_broad.hubness_penalty = "none"
    args_broad.top_k = 7
    args_broad.gmt_topk_list = "7"
    args_broad.gmt_min_genes = 1
    args_broad.emit_small_gene_sets = True
    morphology_profile_query.run(args_broad)
    genes_broad = _genes(Path(args_broad.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv")

    assert "NTRK1" in genes_default
    if "NTRK1" in genes_broad:
        assert genes_default.index("NTRK1") < genes_broad.index("NTRK1")
    assert genes_default[0] == "NTRK1"


def test_morphology_specificity_benchmark_gene_recurrence_penalty_reduces_generic_recurrent_gene(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    refs.write_text(
        "perturbation_id\tf1\tf2\n"
        "A\t1.0\t0.0\n"
        "B\t0.99\t0.01\n"
        "C\t0.98\t0.02\n"
        "SPEC\t0.96\t0.04\n",
        encoding="utf-8",
    )
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\thub_score\tqc_weight\tis_control\n"
        "A\tcompound\tA\t0.1\t1.0\tfalse\n"
        "B\tcompound\tB\t0.1\t1.0\tfalse\n"
        "C\tcompound\tC\t0.1\t1.0\tfalse\n"
        "SPEC\tcompound\tSPEC\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "compound_id\tgene_symbol\tweight\n"
        "A\tADA\t1.0\n"
        "B\tADA\t1.0\n"
        "C\tADA\t1.0\n"
        "SPEC\tKCNN4\t1.0\n",
        encoding="utf-8",
    )
    args_pen = Args()
    args_pen.query_profiles_tsv = str(query)
    args_pen.query_metadata_tsv = None
    args_pen.group_query_by = None
    args_pen.reference_profiles_tsv = str(refs)
    args_pen.reference_metadata_tsv = str(ref_meta)
    args_pen.compound_targets_tsv = str(targets)
    args_pen.feature_stats_tsv = None
    args_pen.feature_schema_tsv = None
    args_pen.out_dir = str(tmp_path / "pen")
    args_pen.polarity = "similar"
    args_pen.top_k = 2
    args_pen.gmt_topk_list = "2"
    args_pen.gmt_min_genes = 1
    args_pen.emit_small_gene_sets = True
    morphology_profile_query.run(args_pen)
    genes_pen_path = Path(args_pen.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv"
    genes_pen = _genes(genes_pen_path)
    scores_pen = _gene_scores(genes_pen_path)

    args_none = Args()
    args_none.query_profiles_tsv = str(query)
    args_none.query_metadata_tsv = None
    args_none.group_query_by = None
    args_none.reference_profiles_tsv = str(refs)
    args_none.reference_metadata_tsv = str(ref_meta)
    args_none.compound_targets_tsv = str(targets)
    args_none.feature_stats_tsv = None
    args_none.feature_schema_tsv = None
    args_none.out_dir = str(tmp_path / "none")
    args_none.polarity = "similar"
    args_none.gene_recurrence_penalty = "none"
    args_none.top_k = 2
    args_none.gmt_topk_list = "2"
    args_none.gmt_min_genes = 1
    args_none.emit_small_gene_sets = True
    morphology_profile_query.run(args_none)
    genes_none_path = Path(args_none.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv"
    genes_none = _genes(genes_none_path)
    scores_none = _gene_scores(genes_none_path)

    assert genes_none[0] == "ADA"
    assert "KCNN4" in genes_pen
    assert scores_pen["KCNN4"] > scores_none["KCNN4"]
    assert scores_pen["ADA"] < scores_none["ADA"]


def test_morphology_direct_target_mode_rescues_distributed_same_target_support(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    query_meta = tmp_path / "query_meta.tsv"
    query_meta.write_text("sample_id\tgene_symbol\nQ1\tNTRK1\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    ref_lines = ["perturbation_id\tf1\tf2"]
    meta_lines = ["perturbation_id\tperturbation_type\tcompound_id\thub_score\tqc_weight\tis_control"]
    target_lines = ["compound_id\tgene_symbol\tweight"]
    ref_lines.append("OFF1\t1.0\t0.0")
    meta_lines.append("OFF1\tcompound\tOFF1\t0.1\t1.0\tfalse")
    target_lines.append("OFF1\tOFF_TARGET\t1.0")
    for idx in range(1, 6):
        ref_lines.append(f"T{idx}\t0.72\t0.28")
        meta_lines.append(f"T{idx}\tcompound\tT{idx}\t0.1\t1.0\tfalse")
        target_lines.append(f"T{idx}\tNTRK1\t1.0")
    refs.write_text("\n".join(ref_lines) + "\n", encoding="utf-8")
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text("\n".join(meta_lines) + "\n", encoding="utf-8")
    targets = tmp_path / "targets.tsv"
    targets.write_text("\n".join(target_lines) + "\n", encoding="utf-8")

    args = Args()
    args.query_profiles_tsv = str(query)
    args.query_metadata_tsv = str(query_meta)
    args.group_query_by = None
    args.reference_profiles_tsv = str(refs)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "direct_target")
    args.mode = "direct_target"
    args.polarity = "similar"
    args.top_k = 2
    args.gmt_topk_list = "2"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    genes = _genes(Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv")
    assert genes[0] == "NTRK1"


def test_morphology_mechanism_mode_backs_off_to_family_support(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    refs.write_text(
        "perturbation_id\tf1\tf2\n"
        "GENERIC\t1.0\t0.0\n"
        "FAM1\t0.78\t0.22\n"
        "FAM2\t0.76\t0.24\n",
        encoding="utf-8",
    )
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\thub_score\tqc_weight\tis_control\n"
        "GENERIC\tcompound\tGENERIC\t0.1\t1.0\tfalse\n"
        "FAM1\tcompound\tFAM1\t0.1\t1.0\tfalse\n"
        "FAM2\tcompound\tFAM2\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "compound_id\tgene_symbol\tweight\n"
        "GENERIC\tOPRM1\t1.0\n"
        "FAM1\tKCNN2\t1.0\n"
        "FAM2\tKCNMA1\t1.0\n",
        encoding="utf-8",
    )
    annotations = tmp_path / "target_annotations.tsv"
    annotations.write_text(
        "gene_symbol\ttarget_family\ttarget_class\tmechanism_label\tpathway_seed\n"
        "OPRM1\tGPCR\tGPCR\tNeuroactive receptor\tN0\n"
        "KCNN2\tIon channel\tPotassium channel\tChannel signaling\tK1\n"
        "KCNMA1\tIon channel\tPotassium channel\tChannel signaling\tK2\n"
        "KCNN4\tIon channel\tPotassium channel\tChannel signaling\tK3\n"
        "KCNK1\tIon channel\tPotassium channel\tChannel signaling\tK4\n",
        encoding="utf-8",
    )

    args = Args()
    args.query_profiles_tsv = str(query)
    args.query_metadata_tsv = None
    args.group_query_by = None
    args.reference_profiles_tsv = str(refs)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.target_annotations_tsv = str(annotations)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "mechanism")
    args.mode = "mechanism"
    args.polarity = "similar"
    args.top_k = 4
    args.gmt_topk_list = "4"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    genes = _genes(Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv")
    assert genes[0] in {"KCNN2", "KCNMA1", "KCNN4", "KCNK1"}
    assert "OPRM1" not in genes[:2]


def test_morphology_mechanism_mode_rejects_unstable_wrong_family_expansion(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    query_meta = tmp_path / "query_meta.tsv"
    query_meta.write_text("sample_id\tperturbation_type\tgene_symbol\nQ1\torf\tNTRK1\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    refs.write_text(
        "perturbation_id\tf1\tf2\n"
        "ION1\t1.00\t0.00\n"
        "RTK1\t0.88\t0.12\n"
        "RTK2\t0.87\t0.13\n",
        encoding="utf-8",
    )
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\tgene_symbol\thub_score\tqc_weight\tis_control\n"
        "ION1\tcompound\tION1\t\t0.1\t1.0\tfalse\n"
        "RTK1\torf\t\tNTRK1\t0.1\t1.0\tfalse\n"
        "RTK2\torf\t\tRET\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text("compound_id\tgene_symbol\tweight\nION1\tKCNMA1\t1.0\n", encoding="utf-8")
    annotations = tmp_path / "target_annotations.tsv"
    annotations.write_text(
        "gene_symbol\ttarget_family\ttarget_class\tmechanism_label\tpathway_seed\n"
        "KCNMA1\tIon channel\tPotassium channel\tChannel signaling\tK_channel\n"
        "NTRK1\tKinase\tReceptor tyrosine kinase\tRTK signaling\tRTK\n"
        "RET\tKinase\tReceptor tyrosine kinase\tRTK signaling\tRTK\n",
        encoding="utf-8",
    )

    args = Args()
    args.query_profiles_tsv = str(query)
    args.query_metadata_tsv = str(query_meta)
    args.group_query_by = None
    args.reference_profiles_tsv = str(refs)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.target_annotations_tsv = str(annotations)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "mechanism_reject_wrong_family")
    args.mode = "mechanism"
    args.polarity = "similar"
    args.top_k = 3
    args.gmt_topk_list = "3"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    genes = _genes(Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv")
    meta = Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.meta.json"
    payload = json.loads(meta.read_text(encoding="utf-8"))
    assert genes[0] in {"NTRK1", "RET"}
    assert payload["summary"]["expansion_decision"]["reason"] != "family_unstable_between_raw_and_penalized"
    assert payload["summary"]["top_family"] in {"Ion channel", "Kinase", None}


def test_morphology_hybrid_expands_when_family_support_is_coherent(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    query_meta = tmp_path / "query_meta.tsv"
    query_meta.write_text("sample_id\tperturbation_type\tgene_symbol\nQ1\torf\tKCNN4\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    refs.write_text(
        "perturbation_id\tf1\tf2\n"
        "ORF1\t0.92\t0.08\n"
        "ORF2\t0.90\t0.10\n"
        "CMP1\t0.86\t0.14\n",
        encoding="utf-8",
    )
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\tgene_symbol\thub_score\tqc_weight\tis_control\n"
        "ORF1\torf\t\tKCNN1\t0.1\t1.0\tfalse\n"
        "ORF2\torf\t\tKCNMA1\t0.1\t1.0\tfalse\n"
        "CMP1\tcompound\tCMP1\t\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text("compound_id\tgene_symbol\tweight\nCMP1\tKCNN4\t1.0\n", encoding="utf-8")
    annotations = tmp_path / "target_annotations.tsv"
    annotations.write_text(
        "gene_symbol\ttarget_family\ttarget_class\tmechanism_label\tpathway_seed\n"
        "KCNN1\tIon channel\tPotassium channel\tChannel signaling\tK_channel\n"
        "KCNMA1\tIon channel\tPotassium channel\tChannel signaling\tK_channel\n"
        "KCNN4\tIon channel\tPotassium channel\tChannel signaling\tK_channel\n"
        "KCNK1\tIon channel\tPotassium channel\tChannel signaling\tK_channel\n",
        encoding="utf-8",
    )
    args = Args()
    args.query_profiles_tsv = str(query)
    args.query_metadata_tsv = str(query_meta)
    args.group_query_by = None
    args.reference_profiles_tsv = str(refs)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.target_annotations_tsv = str(annotations)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "hybrid_expand_orf")
    args.mode = "hybrid"
    args.polarity = "similar"
    args.top_k = 4
    args.gmt_topk_list = "4"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    core_genes = _genes(Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.core.tsv")
    expanded_genes = _genes(Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.expanded.tsv")
    meta = Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.meta.json"
    payload = json.loads(meta.read_text(encoding="utf-8"))
    assert len(expanded_genes) > len(core_genes)
    assert payload["summary"]["expansion_decision"]["allow_expansion"] is True
    assert payload["summary"]["expansion_decision"]["chosen_level"] in {"pathway_seed", "target_class", "mechanism_label", "target_family"}


def test_morphology_compound_query_can_expand_from_stable_genetic_family_support(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    query_meta = tmp_path / "query_meta.tsv"
    query_meta.write_text("sample_id\tperturbation_type\tgene_symbol\nQ1\tcompound\tNTRK1\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    refs.write_text(
        "perturbation_id\tf1\tf2\n"
        "ORF1\t0.95\t0.05\n"
        "ORF2\t0.94\t0.06\n"
        "CMP1\t0.90\t0.10\n",
        encoding="utf-8",
    )
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\tgene_symbol\thub_score\tqc_weight\tis_control\n"
        "ORF1\torf\t\tNTRK1\t0.1\t1.0\tfalse\n"
        "ORF2\torf\t\tRET\t0.1\t1.0\tfalse\n"
        "CMP1\tcompound\tCMP1\t\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text("compound_id\tgene_symbol\tweight\nCMP1\tIGF1R\t1.0\n", encoding="utf-8")
    annotations = tmp_path / "target_annotations.tsv"
    annotations.write_text(
        "gene_symbol\ttarget_family\ttarget_class\tmechanism_label\tpathway_seed\n"
        "NTRK1\tKinase\tReceptor tyrosine kinase\tRTK signaling\tRTK\n"
        "RET\tKinase\tReceptor tyrosine kinase\tRTK signaling\tRTK\n"
        "IGF1R\tKinase\tReceptor tyrosine kinase\tRTK signaling\tRTK\n"
        "KDR\tKinase\tReceptor tyrosine kinase\tRTK signaling\tRTK\n",
        encoding="utf-8",
    )
    args = Args()
    args.query_profiles_tsv = str(query)
    args.query_metadata_tsv = str(query_meta)
    args.group_query_by = None
    args.reference_profiles_tsv = str(refs)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.target_annotations_tsv = str(annotations)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "compound_expand")
    args.mode = "hybrid"
    args.polarity = "similar"
    args.top_k = 4
    args.gmt_topk_list = "4"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    core_genes = _genes(Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.core.tsv")
    expanded_genes = _genes(Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.expanded.tsv")
    meta = json.loads((Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert "NTRK1" in expanded_genes
    assert len(expanded_genes) > len(core_genes)
    assert meta["summary"]["expansion_decision"]["require_same_modality_support"] is False
    assert meta["summary"]["expansion_decision"]["allow_expansion"] is True


def test_morphology_hybrid_expands_from_mechanism_branch_when_core_is_tiny(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    query_meta = tmp_path / "query_meta.tsv"
    query_meta.write_text("sample_id\tperturbation_type\tgene_symbol\nQ1\tcompound\tNTRK1\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    refs.write_text(
        "perturbation_id\tf1\tf2\n"
        "CMP_STRICT\t0.98\t0.02\n"
        "ORF_A\t0.95\t0.05\n"
        "ORF_B\t0.94\t0.06\n"
        "ORF_C\t0.93\t0.07\n",
        encoding="utf-8",
    )
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\tgene_symbol\thub_score\tqc_weight\tis_control\n"
        "CMP_STRICT\tcompound\tCMP_STRICT\t\t0.1\t1.0\tfalse\n"
        "ORF_A\torf\t\tNTRK1\t0.1\t1.0\tfalse\n"
        "ORF_B\torf\t\tRET\t0.1\t1.0\tfalse\n"
        "ORF_C\torf\t\tKDR\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text("compound_id\tgene_symbol\tweight\nCMP_STRICT\tFGFR1\t1.0\n", encoding="utf-8")
    annotations = tmp_path / "target_annotations.tsv"
    annotations.write_text(
        "gene_symbol\ttarget_family\ttarget_class\tmechanism_label\tpathway_seed\n"
        "FGFR1\tKinase\tReceptor tyrosine kinase\tRTK signaling\tFGF_RTK\n"
        "NTRK1\tKinase\tReceptor tyrosine kinase\tRTK signaling\tNTRK_RTK\n"
        "RET\tKinase\tReceptor tyrosine kinase\tRTK signaling\tRET_RTK\n"
        "KDR\tKinase\tReceptor tyrosine kinase\tRTK signaling\tVEGF_RTK\n",
        encoding="utf-8",
    )
    args = Args()
    args.query_profiles_tsv = str(query)
    args.query_metadata_tsv = str(query_meta)
    args.group_query_by = None
    args.reference_profiles_tsv = str(refs)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.target_annotations_tsv = str(annotations)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "hybrid_branch_split")
    args.mode = "hybrid"
    args.top_k = 4
    args.gmt_topk_list = "4"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    program_dir = Path(args.out_dir) / "program=Q1__polarity=similar"
    core_genes = _genes(program_dir / "geneset.core.tsv")
    expanded_genes = _genes(program_dir / "geneset.expanded.tsv")
    meta = json.loads((program_dir / "geneset.meta.json").read_text(encoding="utf-8"))
    assert len(core_genes) <= 2
    assert len(expanded_genes) > len(core_genes)
    assert meta["summary"]["mechanism_branch_neighbor_count"] > meta["summary"]["core_branch_neighbor_count"]
    assert meta["summary"]["expansion_decision"]["allow_expansion"] is True


def test_morphology_retained_label_can_beat_raw_label_mismatch(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    query_meta = tmp_path / "query_meta.tsv"
    query_meta.write_text("sample_id\tperturbation_type\tgene_symbol\nQ1\tcompound\tNTRK1\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    refs.write_text(
        "perturbation_id\tf1\tf2\n"
        "RAW1\t1.00\t0.00\n"
        "RAW2\t0.99\t0.01\n"
        "RTK1\t0.95\t0.05\n"
        "RTK2\t0.94\t0.06\n",
        encoding="utf-8",
    )
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\tgene_symbol\thub_score\tqc_weight\tis_control\n"
        "RAW1\tcompound\tRAW1\t\t9.0\t1.0\tfalse\n"
        "RAW2\tcompound\tRAW2\t\t8.0\t1.0\tfalse\n"
        "RTK1\torf\t\tNTRK1\t0.1\t1.0\tfalse\n"
        "RTK2\torf\t\tRET\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "compound_id\tgene_symbol\tweight\n"
        "RAW1\tOPRM1\t1.0\n"
        "RAW2\tP2RY12\t1.0\n",
        encoding="utf-8",
    )
    annotations = tmp_path / "target_annotations.tsv"
    annotations.write_text(
        "gene_symbol\ttarget_family\ttarget_class\tmechanism_label\tpathway_seed\n"
        "OPRM1\tGPCR\tGPCR\tNeuroactive receptor\tGPCR_seed\n"
        "P2RY12\tGPCR\tGPCR\tNeuroactive receptor\tGPCR_seed\n"
        "NTRK1\tKinase\tReceptor tyrosine kinase\tRTK signaling\tNTRK_RTK\n"
        "RET\tKinase\tReceptor tyrosine kinase\tRTK signaling\tRET_RTK\n",
        encoding="utf-8",
    )
    args = Args()
    args.query_profiles_tsv = str(query)
    args.query_metadata_tsv = str(query_meta)
    args.group_query_by = None
    args.reference_profiles_tsv = str(refs)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.target_annotations_tsv = str(annotations)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "retained_beats_raw")
    args.mode = "mechanism"
    args.top_k = 4
    args.gmt_topk_list = "4"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    meta = json.loads((Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["summary"]["expansion_decision"]["raw_retained_mismatch"] is True
    assert meta["summary"]["expansion_decision"]["allow_expansion"] is True
    assert meta["summary"]["expansion_decision"]["chosen_label"] in {"Receptor tyrosine kinase", "RTK signaling", "Kinase"}


def test_morphology_narrow_label_beats_broad_label(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    refs.write_text(
        "perturbation_id\tf1\tf2\n"
        "A\t0.96\t0.04\n"
        "B\t0.95\t0.05\n",
        encoding="utf-8",
    )
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\thub_score\tqc_weight\tis_control\n"
        "A\tcompound\tA\t0.1\t1.0\tfalse\n"
        "B\tcompound\tB\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text("compound_id\tgene_symbol\tweight\nA\tKCNN1\t1.0\nB\tKCNMA1\t1.0\n", encoding="utf-8")
    annotations = tmp_path / "target_annotations.tsv"
    annotations.write_text(
        "gene_symbol\ttarget_family\ttarget_class\tmechanism_label\tpathway_seed\n"
        "KCNN1\tIon channel\tPotassium channel\tChannel signaling\tPotassium_conductance\n"
        "KCNMA1\tIon channel\tPotassium channel\tChannel signaling\tPotassium_conductance\n"
        "KCNN4\tIon channel\tPotassium channel\tChannel signaling\tPotassium_conductance\n"
        "SCN5A\tIon channel\tSodium channel\tChannel signaling\tSodium_conductance\n",
        encoding="utf-8",
    )
    args = Args()
    args.query_profiles_tsv = str(query)
    args.query_metadata_tsv = None
    args.group_query_by = None
    args.reference_profiles_tsv = str(refs)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.target_annotations_tsv = str(annotations)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "narrow_wins")
    args.mode = "mechanism"
    args.top_k = 4
    args.gmt_topk_list = "4"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    meta = json.loads((Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["summary"]["expansion_decision"]["chosen_level"] in {"pathway_seed", "target_class"}
    expanded_genes = _genes(Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv")
    assert "SCN5A" not in expanded_genes


def test_morphology_expansion_uses_bundle_gene_universe_not_reference_ids(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    refs.write_text("perturbation_id\tf1\tf2\nA\t0.95\t0.05\nB\t0.94\t0.06\nC\t0.10\t0.90\n", encoding="utf-8")
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\thub_score\tqc_weight\tis_control\n"
        "A\tcompound\tA\t0.1\t1.0\tfalse\n"
        "B\tcompound\tB\t0.1\t1.0\tfalse\n"
        "C\tcompound\tC\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text("compound_id\tgene_symbol\tweight\nA\tKCNN1\t1.0\nB\tKCNMA1\t1.0\nC\tKCNN4\t1.0\n", encoding="utf-8")
    annotations = tmp_path / "target_annotations.tsv"
    annotations.write_text(
        "gene_symbol\ttarget_family\ttarget_class\tmechanism_label\tpathway_seed\n"
        "KCNN1\tIon channel\tPotassium channel\tChannel signaling\tPotassium_conductance\n"
        "KCNMA1\tIon channel\tPotassium channel\tChannel signaling\tPotassium_conductance\n"
        "KCNN4\tIon channel\tPotassium channel\tChannel signaling\tPotassium_conductance\n",
        encoding="utf-8",
    )
    args = Args()
    args.query_profiles_tsv = str(query)
    args.query_metadata_tsv = None
    args.group_query_by = None
    args.reference_profiles_tsv = str(refs)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.target_annotations_tsv = str(annotations)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "bundle_gene_universe")
    args.mode = "mechanism"
    args.top_k = 4
    args.gmt_topk_list = "4"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    meta = json.loads((Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.meta.json").read_text(encoding="utf-8"))
    genes = _genes(Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.tsv")
    assert "KCNN4" in genes
    assert meta["summary"]["expansion_decision"]["candidate_scope"] == "bundle"
    assert meta["summary"]["expansion_decision"]["bundle_candidate_genes"] >= 3


def test_morphology_target_nomination_rewards_distributed_support(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    refs.write_text(
        "perturbation_id\tf1\tf2\n"
        "A\t0.98\t0.02\n"
        "B\t0.97\t0.03\n"
        "SHARP\t1.00\t0.00\n",
        encoding="utf-8",
    )
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\thub_score\tqc_weight\tis_control\n"
        "A\tcompound\tA\t0.1\t1.0\tfalse\n"
        "B\tcompound\tB\t0.1\t1.0\tfalse\n"
        "SHARP\tcompound\tSHARP\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "compound_id\tgene_symbol\tweight\n"
        "A\tKCNN4\t1.0\n"
        "B\tKCNN4\t1.0\n"
        "SHARP\tOPRM1\t1.0\n",
        encoding="utf-8",
    )
    args = Args()
    args.query_profiles_tsv = str(query)
    args.query_metadata_tsv = None
    args.group_query_by = None
    args.reference_profiles_tsv = str(refs)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "distributed_nomination")
    args.mode = "direct_target"
    args.top_k = 2
    args.gmt_topk_list = "2"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    meta = json.loads((Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["summary"]["top_target_candidates"][0]["gene_symbol"] == "KCNN4"


def test_morphology_summary_exposes_target_class_and_pathway_seed(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    refs.write_text("perturbation_id\tf1\tf2\nA\t0.95\t0.05\nB\t0.94\t0.06\n", encoding="utf-8")
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\thub_score\tqc_weight\tis_control\n"
        "A\tcompound\tA\t0.1\t1.0\tfalse\n"
        "B\tcompound\tB\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text("compound_id\tgene_symbol\tweight\nA\tKCNN1\t1.0\nB\tKCNMA1\t1.0\n", encoding="utf-8")
    annotations = tmp_path / "target_annotations.tsv"
    annotations.write_text(
        "gene_symbol\ttarget_family\ttarget_class\tmechanism_label\tpathway_seed\n"
        "KCNN1\tIon channel\tPotassium channel\tChannel signaling\tPotassium_conductance\n"
        "KCNMA1\tIon channel\tPotassium channel\tChannel signaling\tPotassium_conductance\n",
        encoding="utf-8",
    )
    args = Args()
    args.query_profiles_tsv = str(query)
    args.query_metadata_tsv = None
    args.group_query_by = None
    args.reference_profiles_tsv = str(refs)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.target_annotations_tsv = str(annotations)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "summary_narrow")
    args.mode = "mechanism"
    args.top_k = 3
    args.gmt_topk_list = "3"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    meta = json.loads((Path(args.out_dir) / "program=Q1__polarity=similar" / "geneset.meta.json").read_text(encoding="utf-8"))
    assert meta["summary"]["top_target_class"] == "Potassium channel"
    assert meta["summary"]["top_pathway_seed"] == "Potassium_conductance"


def test_morphology_hybrid_output_is_distinct_from_mechanism_when_core_exists(tmp_path: Path):
    query = tmp_path / "query.tsv"
    query.write_text("sample_id\tf1\tf2\nQ1\t1.0\t0.0\n", encoding="utf-8")
    query_meta = tmp_path / "query_meta.tsv"
    query_meta.write_text("sample_id\tperturbation_type\tgene_symbol\nQ1\tcompound\tNTRK1\n", encoding="utf-8")
    refs = tmp_path / "refs.tsv"
    refs.write_text(
        "perturbation_id\tf1\tf2\n"
        "CMP1\t0.99\t0.01\n"
        "ORF1\t0.95\t0.05\n"
        "ORF2\t0.94\t0.06\n",
        encoding="utf-8",
    )
    ref_meta = tmp_path / "ref_meta.tsv"
    ref_meta.write_text(
        "perturbation_id\tperturbation_type\tcompound_id\tgene_symbol\thub_score\tqc_weight\tis_control\n"
        "CMP1\tcompound\tCMP1\t\t0.1\t1.0\tfalse\n"
        "ORF1\torf\t\tNTRK1\t0.1\t1.0\tfalse\n"
        "ORF2\torf\t\tRET\t0.1\t1.0\tfalse\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text("compound_id\tgene_symbol\tweight\nCMP1\tFGFR1\t1.0\n", encoding="utf-8")
    annotations = tmp_path / "target_annotations.tsv"
    annotations.write_text(
        "gene_symbol\ttarget_family\ttarget_class\tmechanism_label\tpathway_seed\n"
        "FGFR1\tKinase\tReceptor tyrosine kinase\tRTK signaling\tFGF_RTK\n"
        "NTRK1\tKinase\tReceptor tyrosine kinase\tRTK signaling\tNTRK_RTK\n"
        "RET\tKinase\tReceptor tyrosine kinase\tRTK signaling\tRET_RTK\n",
        encoding="utf-8",
    )
    args = Args()
    args.query_profiles_tsv = str(query)
    args.query_metadata_tsv = str(query_meta)
    args.group_query_by = None
    args.reference_profiles_tsv = str(refs)
    args.reference_metadata_tsv = str(ref_meta)
    args.compound_targets_tsv = str(targets)
    args.target_annotations_tsv = str(annotations)
    args.feature_stats_tsv = None
    args.feature_schema_tsv = None
    args.out_dir = str(tmp_path / "hybrid_merge")
    args.mode = "hybrid"
    args.top_k = 4
    args.gmt_topk_list = "4"
    args.gmt_min_genes = 1
    args.emit_small_gene_sets = True
    morphology_profile_query.run(args)
    program_dir = Path(args.out_dir) / "program=Q1__polarity=similar"
    core_genes = _genes(program_dir / "geneset.core.tsv")
    expanded_genes = _genes(program_dir / "geneset.expanded.tsv")
    merged_genes = _genes(program_dir / "geneset.tsv")
    meta = json.loads((program_dir / "geneset.meta.json").read_text(encoding="utf-8"))
    assert len(core_genes) < len(expanded_genes)
    assert merged_genes != expanded_genes
    assert set(core_genes).issubset(set(merged_genes))
    assert meta["summary"]["hybrid_merge_applied"] is True


def test_morphology_optional_compound_target_confidence_downweights_secondary_target(tmp_path: Path):
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "compound_id\tgene_symbol\tweight\tconfidence\tprimary_target\n"
        "CMP1\tNTRK1\t1.0\t1.0\ttrue\n"
        "CMP1\tOPRM1\t1.0\t0.2\tfalse\n",
        encoding="utf-8",
    )
    from geneset_extractors.extractors.morphology.io import read_compound_targets_auto

    mapping, summary = read_compound_targets_auto(targets)
    assert mapping["CMP1"]["NTRK1"] > mapping["CMP1"]["OPRM1"]
    assert summary["n_rows_using_optional_confidence"] == 2
