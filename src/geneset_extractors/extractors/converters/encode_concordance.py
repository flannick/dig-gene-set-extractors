"""ENCODE accessibility × RNA-seq / GTEx concordance gene sets.

Five library types:
  - ENCODE_accessibility_expression_concordance_raw    (accessible_raw ∩/- rnaseq)
  - ENCODE_accessibility_expression_concordance_control (accessible_ctrl ∩/- rnaseq)
  - ENCODE_background_sneak_delta                      ((raw - ctrl) ∩/- rnaseq)
  - GTEx_x_accessibility_concordance_raw               (gtex_enriched ∩/- accessible_raw)
  - GTEx_x_accessibility_concordance_control           (gtex_enriched ∩/- accessible_ctrl)

Inputs:
  --accessible_genes_dir   Per-biosample raw accessible gene sets (encode_accessible_genes output).
  --contrast_dir           LOO-contrast accessible gene sets (encode_accessibility_contrast output;
                           Up sets are specifically accessible).
  --rnaseq_dir             Per-biosample RNA-seq expressed gene sets.
  --gtex_enriched_dir      Per-tissue GTEx tissue-enriched gene sets (gtex_tissue_enriched output).
  --assay                  Assay label (default: ATAC-seq).

DNase note: DNase accessible gene sets are not currently in V4. To produce DNase concordance sets,
first run encode_accessible_genes with --assay DNase-seq, then encode_accessibility_contrast for
the DNase control sets, then re-run this converter with those outputs.
"""
from __future__ import annotations

import json
from pathlib import Path

from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.core.metadata import make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context

ENCODE_CITATION = (
    "ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the "
    "human genome. Nature. 2012;489(7414):57-74. doi:10.1038/nature11247. "
    "ENCODE/NHGRI; public data. https://www.encodeproject.org"
)
GTEX_CITATION = (
    "GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across "
    "human tissues. Science. 2020;369(6509):1318-1330. doi:10.1126/science.aaz1776. "
    "NIH Common Fund GTEx V8; public aggregate data."
)

# Derived from V3 ENCODE_GTEx_concordance geneset.meta.json biosample/gtex_tissue pairs.
BIOSAMPLE_TO_GTEX_TISSUE: dict[str, str] = {
    "EL": "Cells_EBV_transformed_lymphocytes",
    "adrenal gland": "Adrenal_Gland",
    "amygdala": "Brain_Amygdala",
    "anterior cingulate cortex": "Brain_Anterior_cingulate_cortex__BA24_",
    "aorta": "Artery_Aorta",
    "body of pancreas": "Pancreas",
    "brain": "Brain_Hypothalamus",
    "cerebellum": "Brain_Cerebellum",
    "epithelial cell of prostate": "Prostate",
    "esophagus mucosa": "Esophagus_Mucosa",
    "esophagus muscularis mucosa": "Esophagus_Muscularis",
    "fallopian tube": "Fallopian_Tube",
    "fibroblast of lung": "Lung",
    "frontal cortex": "Brain_Frontal_Cortex__BA9_",
    "heart": "Heart_Atrial_Appendage",
    "heart left ventricle": "Heart_Left_Ventricle",
    "hypothalamus": "Brain_Hypothalamus",
    "kidney": "Kidney_Cortex",
    "left lobe of liver": "Liver",
    "left lung": "Lung",
    "liver": "Liver",
    "lower lobe of left lung": "Lung",
    "lower lobe of right lung": "Lung",
    "lung": "Lung",
    "lung microvascular endothelial cell": "Lung",
    "mucosa of gallbladder": "Bladder",
    "mucosa of urinary bladder": "Bladder",
    "nucleus accumbens": "Brain_Nucleus_accumbens__basal_ganglia_",
    "ovary": "Ovary",
    "pancreas": "Pancreas",
    "progenitor cell of endocrine pancreas": "Pancreas",
    "prostate gland": "Prostate",
    "putamen": "Brain_Putamen__basal_ganglia_",
    "right lobe of liver": "Liver",
    "right lung": "Lung",
    "small intestine": "Small_Intestine_Terminal_Ileum",
    "spinal cord": "Brain_Spinal_cord__cervical_c_1_",
    "spleen": "Spleen",
    "stomach": "Stomach",
    "substantia nigra": "Brain_Substantia_nigra",
    "testis": "Testis",
    "thyroid gland": "Thyroid",
    "upper lobe of left lung": "Lung",
    "upper lobe of right lung": "Lung",
    "urinary bladder": "Bladder",
    "uterus": "Uterus",
    "vagina": "Vagina",
}


def _safe_name(s: str, max_len: int = 80) -> str:
    return "".join(c if c.isalnum() else "_" for c in s)[:max_len]


def _load_genes_from_tsv(path: Path) -> set[str]:
    genes: set[str] = set()
    with path.open("r", encoding="utf-8") as fh:
        next(fh, None)
        for line in fh:
            g = line.strip().split("\t")[0]
            if g:
                genes.add(g)
    return genes


def _parse_cmd_arg(cmd: list[str], flag: str) -> str:
    for i, tok in enumerate(cmd):
        if tok == flag and i + 1 < len(cmd):
            return cmd[i + 1]
    return ""


def _get_biosample(meta: dict) -> str:
    """Extract biosample from V3 top-level or V4 execution command."""
    bs = meta.get("biosample", "")
    if bs:
        return bs
    cmd = meta.get("converter", {}).get("execution", {}).get("command", [])
    return _parse_cmd_arg(cmd, "--biosample")


def _get_direction(meta: dict) -> str:
    """Extract direction from V3 top-level or V4 execution command."""
    d = meta.get("direction", "")
    if d:
        return d
    cmd = meta.get("converter", {}).get("execution", {}).get("command", [])
    return _parse_cmd_arg(cmd, "--direction")


def _load_raw_accessible(src: Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    """Return ({biosample: genes}, {dir_name: biosample}) from accessible-genes dir."""
    by_biosample: dict[str, set[str]] = {}
    dirmap: dict[str, str] = {}
    for sub in sorted(src.iterdir()):
        tsv = sub / "geneset.tsv"
        if not tsv.is_file():
            continue
        genes = _load_genes_from_tsv(tsv)
        if not genes:
            continue
        biosample = sub.name
        meta_path = sub / "geneset.meta.json"
        if meta_path.is_file():
            try:
                bs = _get_biosample(json.loads(meta_path.read_text()))
                if bs:
                    biosample = bs
            except Exception:
                pass
        by_biosample[biosample] = genes
        dirmap[sub.name] = biosample
    return by_biosample, dirmap


def _load_contrast_up(src: Path, dirmap: dict[str, str]) -> dict[str, set[str]]:
    """Return {biosample: genes} for Up (specifically accessible) sets from contrast dir."""
    result: dict[str, set[str]] = {}
    for sub in sorted(src.iterdir()):
        tsv = sub / "geneset.tsv"
        if not tsv.is_file():
            continue
        # Detect direction from meta (V3 or V4) or dir-name suffix
        is_up = sub.name.endswith("_Up")
        biosample_label = sub.name
        meta_path = sub / "geneset.meta.json"
        if meta_path.is_file():
            try:
                meta = json.loads(meta_path.read_text())
                direction = _get_direction(meta)
                if direction:
                    is_up = direction == "Up"
                bs = _get_biosample(meta)
                if bs:
                    biosample_label = bs
            except Exception:
                pass
        if not is_up:
            continue
        genes = _load_genes_from_tsv(tsv)
        if not genes:
            continue
        # biosample_label from V4 meta is the setdir name (e.g. "ENCODE_ATAC_seq_A549_accessible")
        # resolve to actual biosample via the dirmap built from accessible_genes_dir
        biosample = dirmap.get(biosample_label, biosample_label)
        result[biosample] = genes
    return result


def _load_rnaseq(src: Path) -> dict[str, set[str]]:
    result: dict[str, set[str]] = {}
    for sub in sorted(src.iterdir()):
        tsv = sub / "geneset.tsv"
        if not tsv.is_file():
            continue
        genes = _load_genes_from_tsv(tsv)
        if not genes:
            continue
        biosample = sub.name
        meta_path = sub / "geneset.meta.json"
        if meta_path.is_file():
            try:
                bs = _get_biosample(json.loads(meta_path.read_text()))
                if bs:
                    biosample = bs
            except Exception:
                pass
        result[biosample] = genes
    return result


def _load_gtex_enriched(src: Path) -> dict[str, set[str]]:
    result: dict[str, set[str]] = {}
    for sub in sorted(src.iterdir()):
        tsv = sub / "geneset.tsv"
        if not tsv.is_file():
            continue
        genes = _load_genes_from_tsv(tsv)
        if not genes:
            continue
        tissue = sub.name.replace("GTEx_tissue_enriched_", "")
        result[tissue] = genes
    return result


def _write_set(
    out_dir: Path,
    set_name: str,
    genes: set[str],
    parameters: dict,
    description: str,
    assay: str,
) -> None:
    set_dir = out_dir / set_name
    set_dir.mkdir(parents=True, exist_ok=True)
    gene_list = sorted(genes)
    with (set_dir / "geneset.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\n")
        for g in gene_list:
            fh.write(f"{g}\n")
    write_gmt([(set_name, gene_list)], set_dir / "genesets.gmt")
    meta = make_metadata(
        converter_name="encode_concordance",
        parameters=parameters,
        data_type="chromatin_accessibility_x_expression",
        assay=assay,
        organism="human",
        genome_build="GRCh38",
        files=[],
        gene_annotation={"mode": "intersection", "source": "encode_concordance"},
        weights={
            "weight_type": "binary",
            "normalization": {"method": "none", "target_sum": None},
            "aggregation": "set_intersection",
        },
        summary={
            "n_input_features": 0,
            "n_genes": len(gene_list),
            "n_features_assigned": len(gene_list),
            "fraction_features_assigned": 1.0,
            "n_sets_emitted": 1,
        },
        gene_set_description=description,
        output_files=[
            {"path": "genesets.gmt", "role": "gmt_library"},
            {"path": "geneset.tsv", "role": "selected_program"},
            {"path": "geneset.meta.json", "role": "metadata_json"},
        ],
    )
    write_metadata(set_dir / "geneset.meta.json", meta)


def run(args) -> dict[str, object]:
    activate_runtime_context(
        "encode_concordance",
        getattr(args, "provenance_overlay_json", None),
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    assay = getattr(args, "assay", "ATAC-seq")
    safe_assay = _safe_name(assay.replace("-seq", "").replace("-", ""))

    # Load raw accessible genes
    raw_acc: dict[str, set[str]] = {}
    dirmap: dict[str, str] = {}
    accessible_src = ""
    acc_dir_arg = getattr(args, "accessible_genes_dir", None)
    if acc_dir_arg:
        src = Path(acc_dir_arg)
        if not src.is_dir():
            raise NotADirectoryError(f"accessible_genes_dir not found: {src}")
        raw_acc, dirmap = _load_raw_accessible(src)
        accessible_src = str(src)

    # Load control accessible genes (Up sets)
    ctrl_acc: dict[str, set[str]] = {}
    contrast_src = ""
    contrast_dir_arg = getattr(args, "contrast_dir", None)
    if contrast_dir_arg:
        src = Path(contrast_dir_arg)
        if not src.is_dir():
            raise NotADirectoryError(f"contrast_dir not found: {src}")
        ctrl_acc = _load_contrast_up(src, dirmap)
        contrast_src = str(src)

    # Load RNA-seq expressed genes
    rnaseq: dict[str, set[str]] = {}
    rnaseq_src = ""
    rnaseq_dir_arg = getattr(args, "rnaseq_dir", None)
    if rnaseq_dir_arg:
        src = Path(rnaseq_dir_arg)
        if not src.is_dir():
            raise NotADirectoryError(f"rnaseq_dir not found: {src}")
        rnaseq = _load_rnaseq(src)
        rnaseq_src = str(src)

    # Load GTEx tissue-enriched genes
    gtex_enriched: dict[str, set[str]] = {}
    gtex_src = ""
    gtex_dir_arg = getattr(args, "gtex_enriched_dir", None)
    if gtex_dir_arg:
        src = Path(gtex_dir_arg)
        if not src.is_dir():
            raise NotADirectoryError(f"gtex_enriched_dir not found: {src}")
        gtex_enriched = _load_gtex_enriched(src)
        gtex_src = str(src)

    if not raw_acc and not ctrl_acc:
        raise ValueError(
            "No accessible gene data loaded. Provide --accessible_genes_dir and/or --contrast_dir."
        )

    n_sets = 0
    counts: dict[str, int] = {}

    # ── ENCODE × RNA-seq concordance ─────────────────────────────────────────
    if rnaseq:
        for mode, acc_map, acc_src in [
            ("raw", raw_acc, accessible_src),
            ("control", ctrl_acc, contrast_src),
        ]:
            lib = f"ENCODE_accessibility_expression_concordance_{mode}"
            for biosample, acc_genes in acc_map.items():
                expr = rnaseq.get(biosample)
                if expr is None:
                    continue
                safe_bs = _safe_name(biosample)
                base_params = {
                    "assay": assay, "biosample": biosample, "control_mode": mode,
                    "library": lib, "accessible_genes_src": acc_src,
                    "rnaseq_src": rnaseq_src, "encode_citation": ENCODE_CITATION,
                }
                for cls, gene_set, desc_frag in [
                    ("concordant_active", acc_genes & expr,
                     f"accessible ({mode}) AND expressed"),
                    ("open_but_silent",   acc_genes - expr,
                     f"accessible ({mode}) but NOT expressed"),
                    ("expressed_but_closed", expr - acc_genes,
                     f"expressed but NOT accessible ({mode})"),
                ]:
                    if not gene_set:
                        continue
                    name = f"ENCODE_{safe_assay}_{mode}_{cls}_{safe_bs}"
                    _write_set(
                        out_dir, name, gene_set,
                        {**base_params, "class": cls},
                        f"Genes {desc_frag} — ENCODE {assay} accessibility ({mode}) × "
                        f"RNA-seq expression in '{biosample}'. {ENCODE_CITATION}",
                        assay,
                    )
                    n_sets += 1
                    counts[f"rnaseq_{mode}_{cls}"] = counts.get(f"rnaseq_{mode}_{cls}", 0) + 1

    # ── Background sneak delta (raw − control) × RNA-seq ─────────────────────
    if rnaseq and raw_acc and ctrl_acc:
        for biosample in set(raw_acc) & set(ctrl_acc):
            expr = rnaseq.get(biosample)
            if expr is None:
                continue
            delta = raw_acc[biosample] - ctrl_acc[biosample]
            if not delta:
                continue
            safe_bs = _safe_name(biosample)
            base_params = {
                "assay": assay, "biosample": biosample,
                "library": "ENCODE_background_sneak_delta",
                "accessible_genes_src": accessible_src, "contrast_src": contrast_src,
                "rnaseq_src": rnaseq_src, "encode_citation": ENCODE_CITATION,
            }
            for cls, gene_set, desc_frag in [
                ("silent",    delta - expr,
                 "accessible in raw data but absorbed by LOO background AND not expressed"),
                ("expressed", delta & expr,
                 "accessible in raw data but absorbed by LOO background AND expressed"),
            ]:
                if not gene_set:
                    continue
                name = f"ENCODE_{safe_assay}_bgsneak_{cls}_{safe_bs}"
                _write_set(
                    out_dir, name, gene_set,
                    {**base_params, "class": f"background_sneak_{cls}"},
                    f"Background-sneak {cls} genes — ENCODE {assay} in '{biosample}': "
                    f"{desc_frag}. {ENCODE_CITATION}",
                    assay,
                )
                n_sets += 1
                counts[f"bgsneak_{cls}"] = counts.get(f"bgsneak_{cls}", 0) + 1

    # ── GTEx × accessibility concordance ─────────────────────────────────────
    if gtex_enriched:
        for biosample, gtex_tissue in BIOSAMPLE_TO_GTEX_TISSUE.items():
            tissue_genes = gtex_enriched.get(gtex_tissue)
            if tissue_genes is None:
                continue
            safe_bs = _safe_name(biosample)
            for mode, acc_map, acc_src in [
                ("raw", raw_acc, accessible_src),
                ("control", ctrl_acc, contrast_src),
            ]:
                acc_genes = acc_map.get(biosample)
                if not acc_genes:
                    continue
                lib = f"GTEx_x_accessibility_concordance_{mode}"
                base_params = {
                    "assay": assay, "biosample": biosample, "gtex_tissue": gtex_tissue,
                    "control_mode": mode, "library": lib,
                    "accessible_genes_src": acc_src, "gtex_enriched_src": gtex_src,
                    "encode_citation": ENCODE_CITATION, "gtex_citation": GTEX_CITATION,
                }
                for cls, gene_set, desc_frag in [
                    ("enriched_not_accessible", tissue_genes - acc_genes,
                     f"GTEx-enriched in {gtex_tissue} but NOT accessible ({mode})"),
                    ("open_not_enriched",       acc_genes - tissue_genes,
                     f"accessible ({mode}) but NOT GTEx-enriched in {gtex_tissue}"),
                ]:
                    if not gene_set:
                        continue
                    name = f"GTEx_{safe_assay}_{mode}_{cls}_{safe_bs}"
                    _write_set(
                        out_dir, name, gene_set,
                        {**base_params, "class": cls},
                        f"Genes {desc_frag} in '{biosample}'. "
                        f"{ENCODE_CITATION}  {GTEX_CITATION}",
                        assay,
                    )
                    n_sets += 1
                    counts[f"gtex_{mode}_{cls}"] = counts.get(f"gtex_{mode}_{cls}", 0) + 1

    return {"n_sets": n_sets, "counts": counts, "out_dir": str(out_dir)}
