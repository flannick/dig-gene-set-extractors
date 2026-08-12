from __future__ import annotations

import collections
import csv
import zipfile
from pathlib import Path

from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context

ENCODE_CITATION = (
    "ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the "
    "human genome. Nature. 2012;489(7414):57-74. doi:10.1038/nature11247. "
    "ENCODE/NHGRI; public data. https://www.encodeproject.org"
)
GTEX_CITATION = (
    "GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across "
    "human tissues. Science. 2020;369(6509):1318-1330. doi:10.1126/science.aaz1776. "
    "NIH Common Fund GTEx V8; public aggregate data. "
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
)


def _safe_name(s: str, max_len: int = 60) -> str:
    return "".join(c if c.isalnum() else "_" for c in s)[:max_len]


def _load_genes_from_tsv(path: Path) -> set[str]:
    genes: set[str] = set()
    with path.open("r", encoding="utf-8") as fh:
        next(fh, None)  # skip header
        for line in fh:
            g = line.strip().split("\t")[0]
            if g:
                genes.add(g)
    return genes


def _load_accessible_from_dir(src: Path) -> tuple[dict[str, set[str]], list[str]]:
    """Return (label -> genes, list_of_sources) from per-biosample set subdirs."""
    per_biosample: dict[str, set[str]] = {}
    sources: list[str] = []
    for sub in sorted(src.iterdir()):
        tsv = sub / "geneset.tsv"
        if not tsv.is_file():
            continue
        genes = _load_genes_from_tsv(tsv)
        if not genes:
            continue
        label = sub.name
        meta_path = sub / "geneset.meta.json"
        if meta_path.is_file():
            try:
                import json
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                label = meta.get("biosample", label) or label
            except Exception:
                pass
        per_biosample[label] = genes
        sources.append(f"{sub} (per-biosample accessible genes)")
    return per_biosample, sources


def _load_accessible_from_zip(zpath: Path) -> tuple[dict[str, set[str]], list[str]]:
    """Return (label -> genes, list_of_sources) from a zip of per-biosample set dirs."""
    import json
    per_biosample: dict[str, set[str]] = {}
    sources: list[str] = []
    with zipfile.ZipFile(zpath, "r") as zf:
        names = zf.namelist()
        for n in names:
            if not n.endswith("geneset.tsv"):
                continue
            raw = zf.read(n).decode("utf-8", "replace").splitlines()
            genes: set[str] = {ln.strip() for ln in raw[1:] if ln.strip()}
            if not genes:
                continue
            setdir = n.split("/")[-2]
            label = setdir
            meta_name = n.rsplit("/", 1)[0] + "/geneset.meta.json"
            if meta_name in names:
                try:
                    meta = json.loads(zf.read(meta_name).decode("utf-8", "replace"))
                    label = meta.get("biosample", label) or label
                except Exception:
                    pass
            per_biosample[label] = genes
            sources.append(f"{zpath}!{setdir} (per-biosample accessible genes)")
    return per_biosample, sources


def _load_gtex_from_enriched_dir(
    enriched_dir: Path,
) -> tuple[dict[str, set[str]], str]:
    """Load tissue -> genes from a directory of gtex_tissue_enriched output subdirs."""
    tissues: dict[str, set[str]] = {}
    for sub in sorted(enriched_dir.iterdir()):
        gp = sub / "geneset.tsv"
        if not gp.is_file():
            continue
        tissue = sub.name
        for prefix in ("GTEx_tissue_enriched_",):
            tissue = tissue.replace(prefix, "")
        tissues[tissue] = _load_genes_from_tsv(gp)
    return tissues, str(enriched_dir)


def _load_gtex_from_tstat(
    tstat_path: Path,
    threshold: float,
) -> tuple[dict[str, set[str]], str]:
    """Load tissue -> genes from a GTEx t-stat TSV, filtering to t >= threshold."""
    tissues: dict[str, set[str]] = {}
    with tstat_path.open("r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        tissue_names = header[1:]
        for _ in range(len(tissue_names)):
            tissues[tissue_names[_]] = set()
        for row in reader:
            if not row or not row[0]:
                continue
            gene = row[0]
            for ti, val_str in enumerate(row[1:]):
                if not val_str or val_str == "NA":
                    continue
                try:
                    if float(val_str) >= threshold:
                        tissues[tissue_names[ti]].add(gene)
                except ValueError:
                    pass
    # Drop empty tissues
    tissues = {t: g for t, g in tissues.items() if g}
    return tissues, str(tstat_path)


def _write_concordance_set(
    set_dir: Path,
    set_name: str,
    description: str,
    genes: list[str],
    assay: str,
    tissue: str,
    consensus_fraction: float,
    n_accessible_biosamples: int,
    n_tissue_enriched: int,
    accessible_genes_dir: str,
    gtex_enriched_dir: str,
) -> None:
    set_dir.mkdir(parents=True, exist_ok=True)
    with (set_dir / "geneset.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\n")
        for g in genes:
            fh.write(f"{g}\n")
    write_gmt([(set_name, genes)], set_dir / "genesets.gmt")

    meta = make_metadata(
        converter_name="encode_gtex_concordance",
        parameters={
            "assay": assay,
            "tissue": tissue,
            "consensus_fraction": consensus_fraction,
            "n_accessible_biosamples": n_accessible_biosamples,
            "n_tissue_enriched_genes": n_tissue_enriched,
            "accessible_genes_dir": accessible_genes_dir,
            "gtex_enriched_dir": gtex_enriched_dir,
            "encode_citation": ENCODE_CITATION,
            "gtex_citation": GTEX_CITATION,
        },
        data_type="chromatin_accessibility_x_expression",
        assay=assay,
        organism="human",
        genome_build="GRCh38",
        files=[],
        gene_annotation={
            "mode": "inherited",
            "source": "encode_accessible_genes_output_and_gtex_tstat",
        },
        weights={
            "weight_type": "binary",
            "normalization": {"method": "none", "target_sum": None},
            "aggregation": "intersection",
        },
        summary={
            "n_input_features": n_accessible_biosamples,
            "n_genes": len(genes),
            "n_features_assigned": len(genes),
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
        "encode_gtex_concordance",
        getattr(args, "provenance_overlay_json", None),
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    assay = getattr(args, "assay", "ATAC-seq")
    consensus_fraction = float(getattr(args, "consensus_fraction", 0.25))
    tstat_threshold = float(getattr(args, "tstat_threshold", 4.0))

    # Load per-biosample accessible genes
    accessible_zip = getattr(args, "accessible_genes_zip", None)
    accessible_dir = getattr(args, "accessible_genes_dir", None)

    if accessible_zip:
        src_path = Path(accessible_zip)
        if not src_path.is_file():
            raise FileNotFoundError(f"accessible_genes_zip not found: {src_path}")
        per_biosample, _ = _load_accessible_from_zip(src_path)
    elif accessible_dir:
        src_path = Path(accessible_dir)
        if not src_path.is_dir():
            raise NotADirectoryError(f"accessible_genes_dir not found: {src_path}")
        per_biosample, _ = _load_accessible_from_dir(src_path)
    else:
        raise ValueError("Provide --accessible_genes_dir or --accessible_genes_zip")

    if not per_biosample:
        raise ValueError("No accessible-gene records loaded from input")

    # Build consensus: gene accessible in >= consensus_fraction of all biosamples
    n = len(per_biosample)
    need = max(2, round(consensus_fraction * n))
    freq: collections.Counter[str] = collections.Counter()
    for genes in per_biosample.values():
        for g in genes:
            freq[g] += 1
    consensus: set[str] = {g for g, c in freq.items() if c >= need}

    # Load GTEx tissue-enriched genes
    gtex_enriched_dir = getattr(args, "gtex_enriched_dir", None)
    gtex_tstat_tsv = getattr(args, "gtex_tstat_tsv", None)

    if gtex_enriched_dir:
        enriched_path = Path(gtex_enriched_dir)
        if not enriched_path.is_dir():
            raise NotADirectoryError(f"gtex_enriched_dir not found: {enriched_path}")
        tissues, gtex_src = _load_gtex_from_enriched_dir(enriched_path)
    elif gtex_tstat_tsv:
        tstat_path = Path(gtex_tstat_tsv)
        if not tstat_path.is_file():
            raise FileNotFoundError(f"gtex_tstat_tsv not found: {tstat_path}")
        tissues, gtex_src = _load_gtex_from_tstat(tstat_path, tstat_threshold)
        gtex_file_rec = input_file_record(str(tstat_path), "gtex_tstat_tsv")
    else:
        raise ValueError("Provide --gtex_enriched_dir or --gtex_tstat_tsv")

    safe_assay = _safe_name(assay)
    n_sets = 0
    gene_counts: list[int] = []

    for tissue, enriched in sorted(tissues.items()):
        inter = sorted(consensus & enriched)
        if not inter:
            continue
        set_name = f"ENCODE_{safe_assay}_consensus_x_GTEx_enriched_{_safe_name(tissue)}"
        description = (
            f"Genes accessible in >= {consensus_fraction} of ENCODE {assay} biosamples "
            f"({n} total, need >= {need}) AND GTEx tissue-enriched in {tissue} "
            f"(t-stat >= {tstat_threshold}). "
            f"Chromatin accessibility consensus intersected with tissue expression specificity. "
            f"GRCh38; ENCODE/NHGRI + GTEx V8/NIH Common Fund; public. {ENCODE_CITATION}. {GTEX_CITATION}"
        )
        _write_concordance_set(
            set_dir=out_dir / set_name,
            set_name=set_name,
            description=description,
            genes=inter,
            assay=assay,
            tissue=tissue,
            consensus_fraction=consensus_fraction,
            n_accessible_biosamples=n,
            n_tissue_enriched=len(enriched),
            accessible_genes_dir=str(src_path),
            gtex_enriched_dir=str(enriched_path),
        )
        gene_counts.append(len(inter))
        n_sets += 1

    return {
        "n_sets": n_sets,
        "n_accessible_biosamples": n,
        "n_consensus_genes": len(consensus),
        "n_gtex_tissues": len(tissues),
        "n_genes_min": min(gene_counts) if gene_counts else 0,
        "n_genes_max": max(gene_counts) if gene_counts else 0,
        "out_dir": str(out_dir),
    }
