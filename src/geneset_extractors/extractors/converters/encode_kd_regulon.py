from __future__ import annotations

import csv
import gzip
import json
import math
import collections
import socket
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import NamedTuple

from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context

socket.setdefaulttimeout(120)

ENCODE_CITATION = (
    "ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the "
    "human genome. Nature. 2012;489(7414):57-74. doi:10.1038/nature11247. "
    "ENCODE/NHGRI; public data. https://www.encodeproject.org"
)
NCBI_GENE_INFO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/"
    "Homo_sapiens.gene_info.gz"
)
ENCODE_BASE = "https://www.encodeproject.org"


class _FileEntry(NamedTuple):
    target: str
    role: str          # "kd" or "ctrl"
    file_accession: str
    href: str
    experiment_accession: str


def _safe_name(s: str, max_len: int = 60) -> str:
    return "".join(c if c.isalnum() else "_" for c in s)[:max_len]


# ── ENSG → symbol ────────────────────────────────────────────────────────────

def _load_ensg_to_symbol(gene_info_gz: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with gzip.open(gene_info_gz, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            symbol = fields[2]
            for xref in fields[5].split("|"):
                if xref.startswith("Ensembl:"):
                    ensg = xref.split(":", 1)[1]
                    mapping[ensg] = symbol
    return mapping


def _maybe_download_gene_info(cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / "Homo_sapiens.gene_info.gz"
    if not dest.exists():
        urllib.request.urlretrieve(NCBI_GENE_INFO_URL, dest)
    return dest


# ── ENCODE API helpers ────────────────────────────────────────────────────────

def _json_get(url: str) -> dict:
    req = urllib.request.Request(
        url, headers={"Accept": "application/json", "User-Agent": "encode-kd-regulon/1.0"}
    )
    with urllib.request.urlopen(req, timeout=600) as resp:
        return json.load(resp)


def _quant_files_for_experiment(exp_id: str) -> list[tuple[str, str]]:
    """Return [(accession, href)] for released GRCh38 gene-quant TSV files."""
    url = f"{ENCODE_BASE}{exp_id}?format=json"
    data = _json_get(url)
    out: list[tuple[str, str]] = []
    for f in data.get("files", []):
        if not isinstance(f, dict):
            continue
        if (
            f.get("output_type", "").startswith("gene quant")
            and f.get("file_format") == "tsv"
            and f.get("assembly") == "GRCh38"
            and f.get("status") == "released"
        ):
            out.append((f["accession"], f["href"]))
    return out


def _query_api(assays: list[str]) -> list[_FileEntry]:
    """Live API query: experiments → per-experiment quant files → FileEntry list."""
    entries: list[_FileEntry] = []
    by_target: dict[str, dict] = collections.defaultdict(
        lambda: {"kd": [], "ctrl": set()}
    )

    for assay in assays:
        params = [
            ("type", "Experiment"),
            ("assay_title", assay),
            ("assembly", "GRCh38"),
            ("status", "released"),
            ("limit", "all"),
            ("format", "json"),
            ("field", "accession"),
            ("field", "@id"),
            ("field", "target.label"),
            ("field", "possible_controls"),
        ]
        url = f"{ENCODE_BASE}/search/?{urllib.parse.urlencode(params)}"
        data = _json_get(url)
        for exp in data.get("@graph", []):
            tgt = exp.get("target")
            if isinstance(tgt, dict):
                tgt = tgt.get("label")
            if not tgt:
                continue
            exp_id = exp["@id"]
            by_target[tgt]["kd"].append(exp_id)
            for ctrl in exp.get("possible_controls", []):
                cid = ctrl.get("@id") if isinstance(ctrl, dict) else ctrl
                if cid:
                    by_target[tgt]["ctrl"].add(cid)

    for tgt, info in sorted(by_target.items()):
        for exp_id in info["kd"]:
            time.sleep(0.2)
            try:
                for acc, href in _quant_files_for_experiment(exp_id):
                    entries.append(
                        _FileEntry(tgt, "kd", acc, href, exp_id)
                    )
            except Exception:
                pass
        for exp_id in sorted(info["ctrl"]):
            time.sleep(0.2)
            try:
                for acc, href in _quant_files_for_experiment(exp_id):
                    entries.append(
                        _FileEntry(tgt, "ctrl", acc, href, exp_id)
                    )
            except Exception:
                pass

    return entries


def _load_manifest(tsv_path: Path) -> list[_FileEntry]:
    """Load FileEntry list from a pre-built manifest TSV.

    Expected columns (tab-separated, with header):
      target  role  file_accession  href  experiment_accession
    """
    entries: list[_FileEntry] = []
    with tsv_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            entries.append(
                _FileEntry(
                    target=row["target"],
                    role=row["role"],
                    file_accession=row["file_accession"],
                    href=row["href"],
                    experiment_accession=row.get("experiment_accession", ""),
                )
            )
    return entries


# ── expression parsing ────────────────────────────────────────────────────────

def _parse_quant_tsv(fp: Path) -> dict[str, float]:
    """Parse a gene-quantification TSV; return {ENSG_noversion: expression}.

    Handles:
    - ENCODE RNA-seq: columns gene_id + TPM (preferred) or FPKM
    - STAR ReadsPerGene: ENSG rows with 4 columns (unstranded/fwd/rev)
    """
    rows: list[list[str]] = []
    hdr: list[str] | None = None
    with fp.open("r", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        first = next(reader, None)
        if first and any(c in first for c in ("gene_id", "TPM", "FPKM")):
            hdr = first
            rows = list(reader)
        else:
            if first:
                rows.append(first)
            rows.extend(reader)

    vals: dict[str, float] = {}
    if hdr:
        gi = hdr.index("gene_id") if "gene_id" in hdr else 0
        col: int | None = None
        if "TPM" in hdr:
            col = hdr.index("TPM")
        elif "FPKM" in hdr:
            col = hdr.index("FPKM")
        if col is None:
            return {}
        for row in rows:
            if len(row) <= col or len(row) <= gi:
                continue
            g = row[gi]
            if not g.startswith("ENSG"):
                continue
            try:
                vals[g.split(".")[0]] = float(row[col])
            except ValueError:
                pass
    else:
        data = [r for r in rows if r and r[0].startswith("ENSG") and len(r) >= 4]
        if not data:
            return {}
        totals = [0, 0, 0]
        for r in data:
            for j in range(3):
                try:
                    totals[j] += int(float(r[j + 1]))
                except ValueError:
                    pass
        col = totals.index(max(totals)) + 1
        lib = max(totals) or 1
        for r in data:
            try:
                vals[r[0].split(".")[0]] = float(r[col]) / lib * 1e6
            except ValueError:
                pass

    return vals


def _mean_expr(
    file_entries: list[tuple[str, str]],
    tmp_dir: Path,
) -> dict[str, float]:
    """Download quant files, average expression, delete files. Returns {} on failure."""
    acc_sums: dict[str, float] = collections.defaultdict(float)
    n = 0
    for acc, href in file_entries:
        fp = tmp_dir / f"{acc}.tsv"
        try:
            urllib.request.urlretrieve(f"{ENCODE_BASE}{href}", fp)
            v = _parse_quant_tsv(fp)
            fp.unlink(missing_ok=True)
            if v:
                n += 1
                for g, x in v.items():
                    acc_sums[g] += x
        except Exception:
            fp.unlink(missing_ok=True)
    if n == 0:
        return {}
    return {g: s / n for g, s in acc_sums.items()}


# ── output writing ────────────────────────────────────────────────────────────

def _write_kd_set(
    set_dir: Path,
    set_name: str,
    description: str,
    genes: list[str],
    assay: str,
    perturbation: str,
    target: str,
    direction: str,
    lfc_threshold: float,
    expr_threshold: float,
    input_file_rec: dict,
) -> None:
    set_dir.mkdir(parents=True, exist_ok=True)
    with (set_dir / "geneset.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\n")
        for g in genes:
            fh.write(f"{g}\n")
    write_gmt([(set_name, genes)], set_dir / "genesets.gmt")

    meta = make_metadata(
        converter_name="encode_kd_regulon",
        parameters={
            "assay": assay,
            "perturbation": perturbation,
            "target": target,
            "direction": direction,
            "lfc_threshold": lfc_threshold,
            "expr_threshold": expr_threshold,
            "ensg_symbol_source": NCBI_GENE_INFO_URL,
            "encode_citation": ENCODE_CITATION,
            "caveat": (
                f"Cell-line {perturbation} differential expression (K562/HepG2; direct + "
                "indirect/secondary effects). NOT tissue-level human biology. "
                "For direct targets, intersect with binding evidence (eCLIP/TF ChIP-seq). "
                "FC-based signature (no formal replicate statistics; no DESeq2/edgeR)."
            ),
        },
        data_type="knockdown_rna_seq",
        assay=assay,
        organism="human",
        genome_build="GRCh38",
        files=[input_file_rec],
        gene_annotation={
            "mode": "ensg_to_symbol",
            "source": NCBI_GENE_INFO_URL,
        },
        weights={
            "weight_type": "binary",
            "normalization": {"method": "none", "target_sum": None},
            "aggregation": "log2fc_threshold",
        },
        summary={
            "n_input_features": 0,
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
        "encode_kd_regulon",
        getattr(args, "provenance_overlay_json", None),
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    lfc = float(getattr(args, "lfc_threshold", 1.0))
    expr_min = float(getattr(args, "expr_threshold", 1.0))
    assay_str = getattr(args, "assay", "shRNA RNA-seq")
    perturbation = getattr(args, "perturbation", "shRNA knockdown")
    lib_prefix = getattr(args, "lib_prefix", None) or f"ENCODE_{_safe_name(assay_str)}_regulon"
    assays = [a.strip() for a in assay_str.split("|") if a.strip()]

    tmp_dir = out_dir / "_tmp_quant"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    ref_dir = out_dir / "references"

    # Load manifest or query API
    manifest_tsv = getattr(args, "encode_kd_manifest_tsv", None)
    if manifest_tsv:
        manifest_path = Path(manifest_tsv)
        if not manifest_path.is_file():
            raise FileNotFoundError(f"encode_kd_manifest_tsv not found: {manifest_path}")
        entries = _load_manifest(manifest_path)
        input_file_rec = input_file_record(str(manifest_path), "encode_kd_manifest_tsv")
        api_url = None
    else:
        params = [
            ("type", "Experiment"),
            ("assay_title", assays[0]),
            ("assembly", "GRCh38"),
            ("status", "released"),
            ("limit", "all"),
            ("format", "json"),
        ]
        api_url = f"{ENCODE_BASE}/search/?{urllib.parse.urlencode(params)}"
        entries = _query_api(assays)
        input_file_rec = {
            "path": api_url,
            "role": "encode_api_query",
            "sha256": None,
            "size_bytes": None,
        }

    if not entries:
        raise ValueError("No KD file entries loaded — check manifest or API query")

    # ENSG → symbol mapping
    gene_info_path_arg = getattr(args, "ncbi_gene_info_gz", None)
    if gene_info_path_arg:
        gene_info_path = Path(gene_info_path_arg)
        if not gene_info_path.is_file():
            raise FileNotFoundError(f"ncbi_gene_info_gz not found: {gene_info_path_path}")
    else:
        gene_info_path = _maybe_download_gene_info(ref_dir)
    ensg2sym = _load_ensg_to_symbol(gene_info_path)

    # Group files by target and role
    by_target: dict[str, dict[str, list[tuple[str, str]]]] = collections.defaultdict(
        lambda: {"kd": [], "ctrl": []}
    )
    for e in entries:
        by_target[e.target][e.role].append((e.file_accession, e.href))

    n_done = n_skip = 0
    up_counts: list[int] = []
    down_counts: list[int] = []

    for target, files in sorted(by_target.items()):
        safe_tgt = _safe_name(target)
        regulated_dir = out_dir / f"{lib_prefix}_{safe_tgt}_regulated"

        # Resume-safe: skip if already written
        if (regulated_dir / "geneset.tsv").is_file():
            n_done += 1
            continue

        kd_files = files["kd"]
        ctrl_files = files["ctrl"]
        if not kd_files or not ctrl_files:
            n_skip += 1
            continue

        kd_expr = _mean_expr(kd_files, tmp_dir)
        ctrl_expr = _mean_expr(ctrl_files, tmp_dir)
        if not kd_expr or not ctrl_expr:
            n_skip += 1
            continue

        up_genes: set[str] = set()
        dn_genes: set[str] = set()
        for ensg in set(kd_expr) | set(ctrl_expr):
            a = kd_expr.get(ensg, 0.0)
            b = ctrl_expr.get(ensg, 0.0)
            if max(a, b) < expr_min:
                continue
            lfc_val = math.log2((a + 1.0) / (b + 1.0))
            sym = ensg2sym.get(ensg)
            if not sym:
                continue
            if lfc_val >= lfc:
                up_genes.add(sym)
            elif lfc_val <= -lfc:
                dn_genes.add(sym)

        reg_genes = up_genes | dn_genes

        cite_suffix = (
            f"ENCODE {perturbation} RNA-seq (KD vs matched control; "
            f"|log2FC| >= {lfc}, max-expr >= {expr_min} TPM); GRCh38; "
            f"ENSG → symbol: NCBI gene_info ({NCBI_GENE_INFO_URL}); "
            f"ENCODE/NHGRI; public. {ENCODE_CITATION}"
        )

        for suffix, genes, direction in [
            ("regulated", sorted(reg_genes), "both"),
            ("up", sorted(up_genes), "up"),
            ("down", sorted(dn_genes), "down"),
        ]:
            if not genes:
                continue
            set_name = f"{lib_prefix}_{safe_tgt}_{suffix}"
            dir_label = {
                "up": "UP (de-repressed / de-repressed on KD)",
                "down": "DOWN (target-dependent / repressed by factor)",
                "both": "differentially expressed (up or down)",
            }[direction]
            _write_kd_set(
                set_dir=out_dir / set_name,
                set_name=set_name,
                description=(
                    f"Genes {dir_label} on {target} {perturbation} "
                    f"(ENCODE; K562/HepG2; {cite_suffix})"
                ),
                genes=genes,
                assay=assay_str,
                perturbation=perturbation,
                target=target,
                direction=direction,
                lfc_threshold=lfc,
                expr_threshold=expr_min,
                input_file_rec=input_file_rec,
            )

        if up_genes:
            up_counts.append(len(up_genes))
        if dn_genes:
            down_counts.append(len(dn_genes))
        n_done += 1

    # Clean up temp dir if empty
    try:
        tmp_dir.rmdir()
    except OSError:
        pass

    all_counts = up_counts + down_counts
    return {
        "n_targets_done": n_done,
        "n_targets_skipped": n_skip,
        "n_up_sets": len(up_counts),
        "n_down_sets": len(down_counts),
        "n_genes_min": min(all_counts) if all_counts else 0,
        "n_genes_max": max(all_counts) if all_counts else 0,
        "out_dir": str(out_dir),
    }
