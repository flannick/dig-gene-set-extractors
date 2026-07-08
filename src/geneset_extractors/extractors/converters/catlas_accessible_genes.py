from __future__ import annotations

import collections
import glob
import gzip
import os
import urllib.request
from pathlib import Path

from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context

# CATlas scATAC-seq — Zhang et al. 2021 Cell (GSE184462; BICCN/NIH)
CATLAS_CITATION = (
    "Zhang K, Hocker JD, Miller M, et al. "
    "A single-cell atlas of chromatin accessibility in the human genome. "
    "Cell. 2021;184(24):5985-6001.e19. doi:10.1016/j.cell.2021.10.024. "
    "GEO accession: GSE184462 (CATlas). NIH Brain Research through Advancing Innovative "
    "Neurotechnologies (BRAIN) Initiative / BICCN; public."
)

REFGENE_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
)

REFGENE_CITATION = (
    "UCSC refGene annotation for GRCh38/hg38. "
    "Haeussler M et al. The UCSC Genome Browser database: 2019 update. "
    "Nucleic Acids Res. 2019;47(D1):D853-D858. doi:10.1093/nar/gky1095. "
    "Downloaded from: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
)


def _safe_name(s: str, max_len: int = 70) -> str:
    return "".join(c if (c.isalnum() or c in "._-") else "_" for c in s)[:max_len]


def _maybe_download(url: str, dest: Path) -> Path:
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = Path(str(dest) + ".part")
        urllib.request.urlretrieve(url, tmp)
        tmp.rename(dest)
    return dest


def _build_promoter_map(
    refgene_gz: Path,
    window: int,
    bin_size: int,
) -> dict[tuple[str, int], set[str]]:
    """Return {(chrom, bin_idx): {gene_symbols}} from UCSC refGene.txt.gz.

    refGene columns (0-based):
      0=bin  1=name  2=chrom  3=strand  4=txStart  5=txEnd  ... 12=name2 (gene symbol)
    TSS = txStart for '+', txEnd for '-' strand.
    """
    bin2genes: dict[tuple[str, int], set[str]] = collections.defaultdict(set)
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
            if not gene:
                continue
            lo = tss - window
            hi = tss + window
            for b in range(lo // bin_size, hi // bin_size + 1):
                bin2genes[(chrom, b)].add(gene)
    return bin2genes


def _covered_bins(bed_gz: Path, bin_size: int) -> set[tuple[str, int]]:
    """Return set of (chrom, bin_idx) covered by any interval in the BED."""
    cov: set[tuple[str, int]] = set()
    with gzip.open(bed_gz, "rt", encoding="utf-8") as fh:
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
                cov.add((chrom, b))
    return cov


def _cell_type_from_filename(path: Path) -> str:
    """Derive a human-readable cell type from an *_Up.bed.gz filename.

    Convention: {cell_type}_Up.bed.gz, with '___' used as a separator
    for names that contain ' / ' (e.g. 'Astrocyte___Oligodendrocyte').
    """
    stem = path.name
    for suffix in ("_Up.bed.gz", "_Up.bed", ".bed.gz", ".bed"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return stem.replace("___", " / ")


def run(args) -> dict[str, object]:
    activate_runtime_context(
        "catlas_accessible_genes",
        getattr(args, "provenance_overlay_json", None),
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    window = int(getattr(args, "promoter_window", 1000))
    bin_size = int(getattr(args, "bin_size", 1000))

    # Resolve refGene: use provided path or download
    if getattr(args, "refgene_gz", None):
        refgene_path = Path(args.refgene_gz)
        if not refgene_path.exists():
            raise FileNotFoundError(f"refgene_gz not found: {refgene_path}")
    else:
        refgene_path = _maybe_download(
            REFGENE_URL,
            out_dir / "references" / "refGene.hg38.txt.gz",
        )

    # Collect BED files
    bed_dir = Path(args.bed_dir)
    if not bed_dir.is_dir():
        raise NotADirectoryError(f"bed_dir is not a directory: {bed_dir}")
    bed_files = sorted(bed_dir.glob("*_Up.bed.gz"))
    if not bed_files:
        # also accept plain .bed.gz without _Up suffix
        bed_files = sorted(bed_dir.glob("*.bed.gz"))
    if not bed_files:
        raise FileNotFoundError(f"No *.bed.gz files found in {bed_dir}")

    # Build promoter map (shared across all cell types)
    bin2genes = _build_promoter_map(refgene_path, window, bin_size)

    refgene_file_rec = input_file_record(str(refgene_path), "refgene_hg38_gz")

    n_sets = 0
    gene_counts: list[int] = []

    for bed_path in bed_files:
        cell_type = _cell_type_from_filename(bed_path)
        safe_ct = _safe_name(cell_type)
        set_name = f"CATLAS_{safe_ct}_accessible_Up"

        cov = _covered_bins(bed_path, bin_size)
        genes = sorted(
            gene
            for bin_key in cov
            if bin_key in bin2genes
            for gene in bin2genes[bin_key]
        )
        # deduplicate (a gene can be hit by multiple bins)
        genes = sorted(set(genes))
        if not genes:
            continue

        description = (
            f"Genes with promoters (TSS±{window}bp, {bin_size}bp bins, GRCh38 refGene) "
            f"overlapping background-prevalence-controlled specifically-accessible regions "
            f"in cell type '{cell_type}' from CATlas human scATAC-seq "
            f"(Zhang et al. Cell 2021; GSE184462; BICCN/NIH; public). "
            f"Input BED = cell-type-specific Up regions after LOO prevalence background "
            f"control (prevalence < 0.25). "
            f"This converter performs the region→gene promoter-overlap step only; "
            f"background control was applied upstream."
        )

        set_dir = out_dir / set_name
        set_dir.mkdir(parents=True, exist_ok=True)

        with (set_dir / "geneset.tsv").open("w", encoding="utf-8", newline="\n") as fh:
            fh.write("gene\n")
            for g in genes:
                fh.write(f"{g}\n")

        write_gmt([(set_name, genes)], set_dir / "genesets.gmt")

        bed_file_rec = input_file_record(str(bed_path), "catlas_accessible_regions_bed_gz")

        meta = make_metadata(
            converter_name="catlas_accessible_genes",
            parameters={
                "cell_type": cell_type,
                "promoter_window_bp": window,
                "bin_size_bp": bin_size,
                "background_control": "LOO_prevalence_lt_0.25 (applied upstream)",
                "refgene_url": REFGENE_URL,
                "genome_assembly": "GRCh38/hg38",
                "catlas_citation": CATLAS_CITATION,
                "refgene_citation": REFGENE_CITATION,
                "input_bed_file": str(bed_path),
                "method": "region_to_gene_promoter_overlap",
                "caveat": (
                    "Promoter accessibility specificity (LOO prevalence background). "
                    "Not a GC/library-normalized read-count differential."
                ),
            },
            data_type="chromatin_accessibility",
            assay="scATAC-seq",
            organism="human",
            genome_build="GRCh38",
            files=[bed_file_rec, refgene_file_rec],
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
                "aggregation": "promoter_region_overlap",
            },
            summary={
                "n_input_features": len(cov),
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
        gene_counts.append(len(genes))
        n_sets += 1

    return {
        "n_groups": n_sets,
        "n_peaks": sum(gene_counts),
        "n_genes": sum(gene_counts),
        "n_genes_min": min(gene_counts) if gene_counts else 0,
        "n_genes_max": max(gene_counts) if gene_counts else 0,
        "out_dir": str(out_dir),
    }
