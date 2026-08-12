from __future__ import annotations

import gzip
import json
import zipfile
import urllib.request
from pathlib import Path

from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context

CLINVAR_VARIANT_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
)
ENCODE_CITATION = (
    "ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the "
    "human genome. Nature. 2012;489(7414):57-74. doi:10.1038/nature11247. "
    "ENCODE/NHGRI; public data. https://www.encodeproject.org"
)
CLINVAR_CITATION = (
    "Landrum MJ et al. ClinVar: improving access to variant interpretations and "
    "supporting evidence. Nucleic Acids Research. 2018;46(D1):D1062-D1067. "
    "doi:10.1093/nar/gkx1153. NCBI/NLM/NIH; public domain. "
    f"{CLINVAR_VARIANT_SUMMARY_URL}"
)

_Record = tuple[str, frozenset[str], str]
# (setdir_name, genes, target_label)


def _safe_name(s: str, max_len: int = 60) -> str:
    return "".join(c if c.isalnum() else "_" for c in s)[:max_len]


# ── ClinVar gene loading ──────────────────────────────────────────────────────

def _derive_clinvar_genes(vs_gz: Path) -> set[str]:
    """Stream ClinVar variant_summary.txt.gz and collect genes with
    Pathogenic/Likely-pathogenic GRCh38 germline variants."""
    genes: set[str] = set()
    with gzip.open(vs_gz, "rt", encoding="utf-8", errors="replace") as fh:
        raw_hdr = fh.readline().rstrip("\n").lstrip("#")
        hdr = raw_hdr.split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        i_gene = ix.get("GeneSymbol")
        i_sig = ix.get("ClinicalSignificance")
        i_asm = ix.get("Assembly")
        if any(v is None for v in (i_gene, i_sig, i_asm)):
            raise ValueError(
                "variant_summary missing expected columns (GeneSymbol, "
                "ClinicalSignificance, Assembly)"
            )
        for line in fh:
            f = line.rstrip("\n").split("\t")
            needed = max(i_gene, i_sig, i_asm) + 1
            if len(f) < needed:
                continue
            if f[i_asm] != "GRCh38":
                continue
            sig = f[i_sig].lower()
            if "pathogenic" not in sig or "conflict" in sig:
                continue
            gene = f[i_gene]
            if not gene or gene in ("-", "") or ";" in gene:
                continue
            genes.add(gene)
    return genes


def _load_clinvar_from_tsv(tsv: Path) -> set[str]:
    """Load pre-computed ClinVar gene set from a geneset.tsv (header: gene)."""
    genes: set[str] = set()
    with tsv.open("r", encoding="utf-8") as fh:
        next(fh, None)
        for line in fh:
            g = line.strip()
            if g:
                genes.add(g)
    return genes


def _maybe_download_variant_summary(cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / "variant_summary.txt.gz"
    if not dest.exists():
        urllib.request.urlretrieve(CLINVAR_VARIANT_SUMMARY_URL, dest)
    return dest


# ── regulon loading ───────────────────────────────────────────────────────────

def _load_records_from_dir(src: Path) -> list[_Record]:
    recs: list[_Record] = []
    for sub in sorted(src.iterdir()):
        tsv = sub / "geneset.tsv"
        if not tsv.is_file():
            continue
        genes: set[str] = set()
        with tsv.open("r", encoding="utf-8") as fh:
            next(fh, None)
            for line in fh:
                g = line.strip()
                if g:
                    genes.add(g)
        if not genes:
            continue
        target = sub.name
        meta_path = sub / "geneset.meta.json"
        if meta_path.is_file():
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                target = meta.get("target", target) or target
            except Exception:
                pass
        recs.append((sub.name, frozenset(genes), target))
    return recs


def _load_records_from_zip(zpath: Path) -> list[_Record]:
    recs: list[_Record] = []
    with zipfile.ZipFile(zpath, "r") as zf:
        names = zf.namelist()
        for n in names:
            if not n.endswith("geneset.tsv"):
                continue
            setdir = n.split("/")[-2]
            raw = zf.read(n).decode("utf-8", "replace").splitlines()
            genes: set[str] = {ln.strip() for ln in raw[1:] if ln.strip()}
            if not genes:
                continue
            target = setdir
            meta_name = n.rsplit("/", 1)[0] + "/geneset.meta.json"
            if meta_name in names:
                try:
                    meta = json.loads(zf.read(meta_name).decode("utf-8", "replace"))
                    target = meta.get("target", target) or target
                except Exception:
                    pass
            recs.append((setdir, frozenset(genes), target))
    return recs


# ── output writing ────────────────────────────────────────────────────────────

def _write_interop_set(
    set_dir: Path,
    set_name: str,
    description: str,
    genes: list[str],
    assay: str,
    target: str | None,
    n_regulon_factors: int,
    n_clinvar_genes: int,
    regulon_dir: str,
    clinvar_file_rec: dict,
) -> None:
    set_dir.mkdir(parents=True, exist_ok=True)
    with (set_dir / "geneset.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\n")
        for g in genes:
            fh.write(f"{g}\n")
    write_gmt([(set_name, genes)], set_dir / "genesets.gmt")

    meta = make_metadata(
        converter_name="interop_regulon_clinvar",
        parameters={
            "assay": assay,
            "target": target,
            "n_regulon_factors": n_regulon_factors,
            "n_clinvar_pathogenic_genes": n_clinvar_genes,
            "regulon_dir": regulon_dir,
            "operation": "regulon_x_clinvar_pathogenic",
            "encode_citation": ENCODE_CITATION,
            "clinvar_citation": CLINVAR_CITATION,
            "note": (
                "Intersection of factor-specific regulon genes with ClinVar pathogenic genes. "
                "Identifies disease genes regulated by or bound to this factor. "
                "Cross-NIH resource (NHGRI ENCODE + NLM ClinVar)."
            ),
        },
        data_type="cross_resource_intersection",
        assay=assay,
        organism="human",
        genome_build="GRCh38",
        files=[clinvar_file_rec],
        gene_annotation={
            "mode": "intersection",
            "source": "regulon_x_clinvar",
        },
        weights={
            "weight_type": "binary",
            "normalization": {"method": "none", "target_sum": None},
            "aggregation": "set_intersection",
        },
        summary={
            "n_input_features": n_regulon_factors,
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


# ── entry point ───────────────────────────────────────────────────────────────

def run(args) -> dict[str, object]:
    activate_runtime_context(
        "interop_regulon_clinvar",
        getattr(args, "provenance_overlay_json", None),
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    assay = getattr(args, "assay", "TF ChIP-seq")
    lib_prefix = getattr(args, "lib_prefix", None) or f"ENCODE_{_safe_name(assay)}_x_ClinVar"
    min_genes = int(getattr(args, "min_genes", 1))
    ref_dir = out_dir / "references"

    # ── Load ClinVar pathogenic genes ─────────────────────────────────────────
    clinvar_tsv = getattr(args, "clinvar_genes_tsv", None)
    variant_summary_gz = getattr(args, "variant_summary_gz", None)

    if clinvar_tsv:
        cv_path = Path(clinvar_tsv)
        if not cv_path.is_file():
            raise FileNotFoundError(f"clinvar_genes_tsv not found: {cv_path}")
        clinvar_genes = _load_clinvar_from_tsv(cv_path)
        clinvar_file_rec = input_file_record(str(cv_path), "clinvar_genes_tsv")
        clinvar_src_label = str(cv_path)
    elif variant_summary_gz:
        vs_path = Path(variant_summary_gz)
        if not vs_path.is_file():
            raise FileNotFoundError(f"variant_summary_gz not found: {vs_path}")
        clinvar_genes = _derive_clinvar_genes(vs_path)
        clinvar_file_rec = input_file_record(str(vs_path), "variant_summary_gz")
        clinvar_src_label = str(vs_path)
    else:
        vs_path = _maybe_download_variant_summary(ref_dir)
        clinvar_genes = _derive_clinvar_genes(vs_path)
        clinvar_file_rec = {
            "path": CLINVAR_VARIANT_SUMMARY_URL,
            "role": "variant_summary_gz",
            "sha256": None,
            "size_bytes": None,
        }
        clinvar_src_label = CLINVAR_VARIANT_SUMMARY_URL

    if not clinvar_genes:
        raise ValueError("ClinVar gene set is empty — check input file or download")

    # ── Load per-factor regulon sets ──────────────────────────────────────────
    regulon_zip = getattr(args, "regulon_zip", None)
    regulon_dir = getattr(args, "regulon_dir", None)

    if regulon_zip:
        rp = Path(regulon_zip)
        if not rp.is_file():
            raise FileNotFoundError(f"regulon_zip not found: {rp}")
        recs = _load_records_from_zip(rp)
        regulon_file_rec = input_file_record(str(rp), "regulon_zip")
    elif regulon_dir:
        rp = Path(regulon_dir)
        if not rp.is_dir():
            raise NotADirectoryError(f"regulon_dir not found: {rp}")
        recs = _load_records_from_dir(rp)
    else:
        raise ValueError("Provide --regulon_dir or --regulon_zip")

    if not recs:
        raise ValueError("No regulon records loaded from input")

    n_factors = len(recs)
    n_clinvar = len(clinvar_genes)

    # ── Per-factor intersection ───────────────────────────────────────────────
    union_disease: set[str] = set()
    n_sets = 0
    gene_counts: list[int] = []

    for _setdir, genes, target in recs:
        inter = sorted(genes & clinvar_genes)
        if len(inter) < min_genes:
            continue
        union_disease.update(inter)
        set_name = f"{lib_prefix}_{_safe_name(target)}_disease_targets"
        _write_interop_set(
            set_dir=out_dir / set_name,
            set_name=set_name,
            description=(
                f"ClinVar Pathogenic/Likely-pathogenic genes specifically targeted "
                f"by {target} ({assay}). "
                f"Intersection of factor regulon ({n_factors} factors total) with "
                f"ClinVar pathogenic gene set ({n_clinvar} genes; GRCh38; "
                f"{clinvar_src_label}). "
                f"Cross-NIH ENCODE x ClinVar interoperability. "
                f"{ENCODE_CITATION}. {CLINVAR_CITATION}"
            ),
            genes=inter,
            assay=assay,
            target=target,
            n_regulon_factors=n_factors,
            n_clinvar_genes=n_clinvar,
            regulon_dir=str(rp),
            clinvar_file_rec=clinvar_file_rec,
        )
        gene_counts.append(len(inter))
        n_sets += 1

    # ── Union: disease genes bound/targeted by any factor ────────────────────
    union_sorted = sorted(union_disease)
    n_union = 0
    if len(union_sorted) >= min_genes:
        union_name = f"{lib_prefix}_any_factor_disease_targets"
        _write_interop_set(
            set_dir=out_dir / union_name,
            set_name=union_name,
            description=(
                f"ClinVar Pathogenic/Likely-pathogenic genes targeted by at least one "
                f"{assay} factor (union across {n_factors} factors). "
                f"ClinVar pathogenic gene set: {n_clinvar} genes; GRCh38; "
                f"{clinvar_src_label}. "
                f"Cross-NIH ENCODE x ClinVar interoperability. "
                f"{ENCODE_CITATION}. {CLINVAR_CITATION}"
            ),
            genes=union_sorted,
            assay=assay,
            target=None,
            n_regulon_factors=n_factors,
            n_clinvar_genes=n_clinvar,
            regulon_dir=str(rp),
            clinvar_file_rec=clinvar_file_rec,
        )
        n_union = 1

    all_counts = gene_counts + ([len(union_sorted)] if n_union else [])
    return {
        "n_regulon_factors": n_factors,
        "n_clinvar_genes": n_clinvar,
        "n_per_factor_sets": n_sets,
        "n_union_sets": n_union,
        "n_union_disease_genes": len(union_disease),
        "n_genes_min": min(gene_counts) if gene_counts else 0,
        "n_genes_max": max(gene_counts) if gene_counts else 0,
        "out_dir": str(out_dir),
    }
