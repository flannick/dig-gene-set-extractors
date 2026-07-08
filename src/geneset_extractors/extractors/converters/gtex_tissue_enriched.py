from __future__ import annotations

import csv
import gzip
import urllib.request
from pathlib import Path

from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context

GTEX_V8_MEDIAN_TPM_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
)


def _safe_name(tissue: str, max_len: int = 60) -> str:
    return "".join(c if c.isalnum() else "_" for c in tissue)[:max_len]


def _maybe_download(path_or_url: str, out_dir: Path) -> Path:
    if path_or_url.startswith(("http://", "https://")):
        dl_dir = out_dir / "downloads"
        dl_dir.mkdir(parents=True, exist_ok=True)
        dest = dl_dir / Path(path_or_url.rstrip("/").split("/")[-1])
        if not dest.exists():
            urllib.request.urlretrieve(path_or_url, dest)
        return dest
    return Path(path_or_url)


def _load_tstat_matrix(path: Path) -> tuple[list[str], list[str], dict[str, list[float]]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header: list[str] = []
        for row in reader:
            if not row or not row[0]:
                continue
            if row[0].startswith("#"):
                continue
            try:
                int(row[0])
                continue  # GCT dimension line
            except ValueError:
                pass
            header = row
            break
        gct = len(header) > 1 and header[1].lower() in ("description", "name", "id")
        val_start = 2 if gct else 1
        gene_col = 1 if gct else 0
        tissues = header[val_start:]
        genes: list[str] = []
        gene_tstat: dict[str, list[float]] = {}
        for row in reader:
            if not row or not row[0]:
                continue
            gene = (row[gene_col] if len(row) > gene_col else row[0]).strip()
            if not gene or gene in ("", "NA"):
                continue
            try:
                values = [float(x) for x in row[val_start:]]
            except ValueError:
                continue
            if len(values) != len(tissues):
                continue
            genes.append(gene)
            gene_tstat[gene] = values
    return tissues, genes, gene_tstat


def run(args) -> dict[str, object]:
    activate_runtime_context("gtex_tissue_enriched", getattr(args, "provenance_overlay_json", None))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tstat_path = _maybe_download(args.gtex_tstat_tsv, out_dir)
    tissues, genes, gene_tstat = _load_tstat_matrix(tstat_path)
    threshold = float(args.tstat_threshold)

    file_rec = input_file_record(str(tstat_path), "gtex_tstat_tsv")

    n_sets = 0
    gene_counts: list[int] = []

    for ti, tissue in enumerate(tissues):
        enriched = sorted(
            [(g, gene_tstat[g][ti]) for g in genes if gene_tstat[g][ti] >= threshold],
            key=lambda kv: -kv[1],
        )
        if not enriched:
            continue
        gene_symbols = [g for g, _ in enriched]
        set_name = f"GTEx_tissue_enriched_{_safe_name(tissue)}"
        description = (
            f"Genes relatively enriched (GTEx t-stat>={threshold}; relative tissue specificity, "
            f"NOT absolute expression) in {tissue}. "
            f"Derived from GTEx V8 median TPM (NIH Common Fund; public aggregate)."
        )

        set_dir = out_dir / set_name
        set_dir.mkdir(parents=True, exist_ok=True)

        with (set_dir / "geneset.tsv").open("w", encoding="utf-8", newline="\n") as fh:
            fh.write("gene\tgtex_tstat\n")
            for g, t in enriched:
                fh.write(f"{g}\t{t:.6f}\n")

        write_gmt([(set_name, gene_symbols)], set_dir / "genesets.gmt")

        meta = make_metadata(
            converter_name="gtex_tissue_enriched",
            parameters={
                "tstat_threshold": threshold,
                "tissue": tissue,
                "gtex_version": "v8",
                "gtex_source_url": GTEX_V8_MEDIAN_TPM_URL,
            },
            data_type="transcriptomics",
            assay="bulk",
            organism="human",
            genome_build="GRCh38",
            files=[file_rec],
            gene_annotation={"mode": "provided", "source": "gtex_hgnc_symbols"},
            weights={
                "weight_type": "nonnegative",
                "normalization": {"method": "none", "target_sum": None},
                "aggregation": "threshold_filter",
            },
            summary={
                "n_input_features": len(genes),
                "n_genes": len(gene_symbols),
                "n_features_assigned": len(gene_symbols),
                "fraction_features_assigned": len(gene_symbols) / len(genes) if genes else 0.0,
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
        gene_counts.append(len(gene_symbols))
        n_sets += 1

    return {
        "n_peaks": sum(gene_counts),
        "n_genes": sum(gene_counts),
        "n_groups": n_sets,
        "n_genes_min": min(gene_counts) if gene_counts else 0,
        "n_genes_max": max(gene_counts) if gene_counts else 0,
        "out_dir": str(out_dir),
    }
