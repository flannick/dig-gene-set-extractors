from __future__ import annotations

import collections
import csv
import gzip
import io
import json
import os
import time
import urllib.parse
import urllib.request
from pathlib import Path

from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context

ENCODE_SEARCH_URL = "https://www.encodeproject.org/search/"
ENCODE_FILE_BASE = "https://www.encodeproject.org"
REFGENE_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
)

REFGENE_CITATION = (
    "UCSC refGene annotation for GRCh38/hg38. "
    "Haeussler M et al. The UCSC Genome Browser database: 2019 update. "
    "Nucleic Acids Res. 2019;47(D1):D853-D858. doi:10.1093/nar/gky1095. "
    f"Downloaded from: {REFGENE_URL}"
)

ENCODE_CITATION = (
    "ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the "
    "human genome. Nature. 2012;489(7414):57-74. doi:10.1038/nature11247. "
    "ENCODE/NHGRI; public data. https://www.encodeproject.org"
)

# Peak output types accepted (in priority order for per-biosample dedup)
_PREFERRED_OUTPUT_TYPES = (
    "IDR thresholded peaks",
    "conservative IDR thresholded peaks",
    "optimal IDR thresholded peaks",
    "pseudoreplicated IDR thresholded peaks",
    "peaks",
)


def _safe_name(s: str, max_len: int = 60) -> str:
    return "".join(c if c.isalnum() else "_" for c in s)[:max_len]


def _maybe_download_refgene(out_dir: Path) -> Path:
    dest = out_dir / "references" / "refGene.hg38.txt.gz"
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = Path(str(dest) + ".part")
        urllib.request.urlretrieve(REFGENE_URL, tmp)
        tmp.rename(dest)
    return dest


def _build_gene_tss_map(
    refgene_gz: Path,
) -> dict[str, set[tuple[str, int]]]:
    """Return {gene_symbol: {(chrom, tss)}} from UCSC refGene.txt.gz."""
    gene_tss: dict[str, set[tuple[str, int]]] = collections.defaultdict(set)
    with gzip.open(refgene_gz, "rt", encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13:
                continue
            chrom, strand = fields[2], fields[3]
            if "_" in chrom or chrom == "chrM":
                continue
            try:
                tss = int(fields[4]) if strand == "+" else int(fields[5])
            except ValueError:
                continue
            gene = fields[12].strip()
            if gene:
                gene_tss[gene].add((chrom, tss))
    return gene_tss


def _accessible_genes(
    peak_gz: Path,
    gene_tss: dict[str, set[tuple[str, int]]],
    window: int,
    bin_size: int,
) -> tuple[list[str], int]:
    """Return (sorted accessible gene list, n_peaks) for one peak BED file."""
    covered: dict[str, set[int]] = collections.defaultdict(set)
    n_peaks = 0
    with gzip.open(peak_gz, "rt", encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            for b in range(start // bin_size, end // bin_size + 1):
                covered[chrom].add(b)
            n_peaks += 1
    acc: list[str] = []
    for gene, tsss in gene_tss.items():
        for chrom, tss in tsss:
            cb = covered.get(chrom)
            if cb and any(
                b in cb
                for b in range((tss - window) // bin_size, (tss + window) // bin_size + 1)
            ):
                acc.append(gene)
                break
    return sorted(set(acc)), n_peaks


def _query_encode_api(
    assay: str,
    output_type: str,
) -> tuple[dict[str, dict], str]:
    """Query ENCODE portal for released peak files.

    Returns (biosample -> {file_accession, experiment_accession, href}, api_url).
    One entry per biosample (first preferred output_type wins).
    """
    params = [
        ("type", "File"),
        ("assay_title", assay),
        ("file_format", "bed"),
        ("assembly", "GRCh38"),
        ("status", "released"),
        ("limit", "all"),
        ("format", "json"),
        ("field", "accession"),
        ("field", "href"),
        ("field", "output_type"),
        ("field", "dataset"),
        ("field", "biosample_ontology"),
    ]
    api_url = ENCODE_SEARCH_URL + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(
        api_url,
        headers={"Accept": "application/json", "User-Agent": "geneset-extractors/1.0"},
    )
    with urllib.request.urlopen(req, timeout=300) as resp:
        data = json.load(resp)

    # Group by biosample, keep best output_type per biosample
    by_biosample: dict[str, dict] = {}
    priority = {ot: i for i, ot in enumerate(_PREFERRED_OUTPUT_TYPES)}

    for f in data.get("@graph", []):
        ot = f.get("output_type", "")
        if ot not in priority:
            continue
        href = f.get("href")
        if not href:
            continue
        bo = f.get("biosample_ontology", {})
        biosample = bo.get("term_name") if isinstance(bo, dict) else None
        if not biosample:
            continue
        exp = f.get("dataset", "").strip("/").split("/")[-1]
        file_acc = f.get("accession", "")
        entry = {
            "file_accession": file_acc,
            "experiment_accession": exp,
            "biosample": biosample,
            "href": href,
            "output_type": ot,
        }
        existing = by_biosample.get(biosample)
        if existing is None or priority[ot] < priority[existing["output_type"]]:
            by_biosample[biosample] = entry

    return by_biosample, api_url


def _load_manifest(path: Path) -> dict[str, dict]:
    """Load a pre-downloaded manifest TSV into the same structure as _query_encode_api."""
    by_biosample: dict[str, dict] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            biosample = row.get("biosample", "").strip()
            if not biosample:
                continue
            by_biosample[biosample] = {
                "file_accession": row.get("file_accession", "").strip(),
                "experiment_accession": row.get("experiment_accession", "").strip(),
                "biosample": biosample,
                "href": row.get("href", "").strip(),
                "output_type": row.get("output_type", "").strip(),
            }
    return by_biosample


def _write_set(
    set_dir: Path,
    set_name: str,
    description: str,
    genes: list[str],
    n_peaks: int,
    assay: str,
    biosample: str,
    file_acc: str,
    exp_acc: str,
    output_type: str,
    api_url: str,
    window: int,
    bin_size: int,
    peak_file_rec: dict,
    refgene_file_rec: dict,
) -> None:
    set_dir.mkdir(parents=True, exist_ok=True)

    with (set_dir / "geneset.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\n")
        for g in genes:
            fh.write(f"{g}\n")

    write_gmt([(set_name, genes)], set_dir / "genesets.gmt")

    cite = (
        f"Derived from ENCODE {assay} peak file {file_acc} "
        f"(experiment {exp_acc}; biosample '{biosample}'; output_type '{output_type}'; "
        f"GRCh38; ENCODE/NHGRI; public). "
        f"Gene assignment: promoter TSS±{window}bp overlap ({bin_size}bp bins), "
        f"UCSC refGene hg38. {ENCODE_CITATION}"
    )

    meta = make_metadata(
        converter_name="encode_accessible_genes",
        parameters={
            "assay": assay,
            "biosample": biosample,
            "encode_file_accession": file_acc,
            "encode_experiment_accession": exp_acc,
            "encode_output_type": output_type,
            "encode_api_query_url": api_url,
            "encode_file_download_url": ENCODE_FILE_BASE + (peak_file_rec.get("url", "")),
            "promoter_window_bp": window,
            "bin_size_bp": bin_size,
            "n_input_peaks": n_peaks,
            "genome_assembly": "GRCh38",
            "refgene_url": REFGENE_URL,
            "refgene_citation": REFGENE_CITATION,
            "encode_citation": ENCODE_CITATION,
            "method": "atac_peak_to_gene_promoter_overlap",
            "caveat": (
                "Single representative peak file per biosample; "
                "not a read-count differential."
            ),
        },
        data_type="chromatin_accessibility",
        assay=assay,
        organism="human",
        genome_build="GRCh38",
        files=[peak_file_rec, refgene_file_rec],
        gene_annotation={
            "mode": "promoter_overlap",
            "source": "UCSC_refGene_hg38",
            "url": REFGENE_URL,
            "window_bp": window,
            "bin_size_bp": bin_size,
            "citation": REFGENE_CITATION,
        },
        weights={
            "weight_type": "nonnegative",
            "normalization": {"method": "none", "target_sum": None},
            "aggregation": "promoter_peak_overlap",
        },
        summary={
            "n_input_features": n_peaks,
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
        "encode_accessible_genes",
        getattr(args, "provenance_overlay_json", None),
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    assay = getattr(args, "assay", "ATAC-seq")
    output_type = getattr(args, "output_type", "IDR thresholded peaks")
    window = int(getattr(args, "promoter_window", 1000))
    bin_size = int(getattr(args, "bin_size", 1000))
    tmp_dir = out_dir / "_peak_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # Resolve refGene
    if getattr(args, "refgene_gz", None):
        refgene_path = Path(args.refgene_gz)
        if not refgene_path.exists():
            raise FileNotFoundError(f"refgene_gz not found: {refgene_path}")
    else:
        refgene_path = _maybe_download_refgene(out_dir)

    refgene_file_rec = input_file_record(str(refgene_path), "refgene_hg38_gz")
    gene_tss = _build_gene_tss_map(refgene_path)

    # Get biosample list: manifest or live API
    manifest_path = getattr(args, "encode_manifest_tsv", None)
    if manifest_path:
        by_biosample = _load_manifest(Path(manifest_path))
        api_url = f"manifest:{manifest_path}"
    else:
        by_biosample, api_url = _query_encode_api(assay, output_type)

    n_done = n_skipped = n_fail = 0
    gene_counts: list[int] = []

    for biosample, entry in sorted(by_biosample.items()):
        href = entry.get("href", "")
        file_acc = entry.get("file_accession", "unknown")
        exp_acc = entry.get("experiment_accession", "unknown")
        ot = entry.get("output_type", output_type)

        safe_bs = _safe_name(biosample)
        set_name = f"ENCODE_{_safe_name(assay)}_{safe_bs}_accessible"
        set_dir = out_dir / set_name

        # Resume-safe: skip if already written
        if (set_dir / "geneset.tsv").exists():
            n_skipped += 1
            continue

        if not href:
            n_fail += 1
            continue

        peak_path = tmp_dir / f"{file_acc}.bed.gz"
        try:
            urllib.request.urlretrieve(ENCODE_FILE_BASE + href, peak_path)
            peak_file_rec = input_file_record(str(peak_path), "encode_peak_bed_gz")
            # Stash the download URL in the record for provenance
            peak_file_rec["url"] = href
            genes, n_peaks = _accessible_genes(peak_path, gene_tss, window, bin_size)
        except Exception as exc:
            n_fail += 1
            try:
                peak_path.unlink(missing_ok=True)
            except OSError:
                pass
            continue
        finally:
            pass

        description = (
            f"Genes with promoter-proximal (TSS±{window}bp, {bin_size}bp bins, GRCh38 refGene) "
            f"ENCODE {assay} accessibility in biosample '{biosample}'. "
            f"Derived from ENCODE peak file {file_acc} (experiment {exp_acc}; "
            f"output_type '{ot}'; GRCh38; ENCODE/NHGRI; public)."
        )

        _write_set(
            set_dir=set_dir,
            set_name=set_name,
            description=description,
            genes=genes,
            n_peaks=n_peaks,
            assay=assay,
            biosample=biosample,
            file_acc=file_acc,
            exp_acc=exp_acc,
            output_type=ot,
            api_url=api_url,
            window=window,
            bin_size=bin_size,
            peak_file_rec=peak_file_rec,
            refgene_file_rec=refgene_file_rec,
        )
        # Delete the downloaded BED after provenance is written
        try:
            peak_path.unlink(missing_ok=True)
        except OSError:
            pass

        gene_counts.append(len(genes))
        n_done += 1

        # Brief pause to be polite to ENCODE portal
        time.sleep(0.5)

    # Clean up tmp dir if empty
    try:
        tmp_dir.rmdir()
    except OSError:
        pass

    return {
        "n_groups": n_done,
        "n_skipped": n_skipped,
        "n_fail": n_fail,
        "n_genes_min": min(gene_counts) if gene_counts else 0,
        "n_genes_max": max(gene_counts) if gene_counts else 0,
        "out_dir": str(out_dir),
    }
