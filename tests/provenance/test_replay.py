import json
from pathlib import Path

from geneset_extractors.cli import main
from geneset_extractors.converters import rna_deg


class Args:
    deg_tsv = "tests/data/toy_deg.tsv"
    out_dir = "tests/tmp/rna_deg"
    organism = "human"
    genome_build = "hg38"
    signature_name = "toy"
    gene_id_column = "gene_id"
    gene_symbol_column = None
    stat_column = None
    logfc_column = None
    padj_column = None
    pvalue_column = None
    score_column = None
    score_mode = "auto"
    padj_max = None
    pvalue_max = None
    min_abs_logfc = None
    duplicate_gene_policy = "max_abs"
    neglog10p_cap = 50.0
    neglog10p_eps = 1e-300
    exclude_gene_regex = None
    disable_default_excludes = False
    gtf = None
    gtf_gene_id_field = "gene_id"
    gtf_source = None
    select = "top_k"
    top_k = 200
    quantile = 0.01
    min_score = 0.0
    normalize = "within_set_l1"
    emit_full = True
    emit_gmt = False
    gmt_out = None
    gmt_prefer_symbol = True
    gmt_require_symbol = False
    gmt_biotype_allowlist = "protein_coding"
    gmt_min_genes = 1
    gmt_max_genes = 10
    gmt_topk_list = "3"
    gmt_mass_list = ""
    gmt_split_signed = True
    gmt_emit_abs = False
    gmt_source = "full"
    emit_small_gene_sets = True
    provenance_overlay_json = None
    emit_replay = True
    emit_compact_provenance = True
    s3_mirror_prefix = None


def test_cli_compact_and_render_replay(tmp_path: Path):
    args = Args()
    args.out_dir = str(tmp_path / "rna_deg_cli")
    rna_deg.run(args)

    prov = Path(args.out_dir) / "geneset.provenance.json"
    compact = tmp_path / "compact.json"
    replay = tmp_path / "reproduce.sh"

    assert main(["provenance", "compact", str(prov), "--out", str(compact)]) == 0
    assert main(["provenance", "render-replay", str(prov), "--out", str(replay)]) == 0

    compact_payload = json.loads(compact.read_text(encoding="utf-8"))
    assert "nodes" in compact_payload
    assert "edges" in compact_payload
    assert "geneset_extractors.cli" in replay.read_text(encoding="utf-8")

