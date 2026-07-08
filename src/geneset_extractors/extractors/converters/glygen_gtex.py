from __future__ import annotations

import csv
import gzip
import io
import urllib.request
from pathlib import Path

from geneset_extractors.core.gmt import write_gmt
from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context

GLYGEN_BASE_URL = "https://data.glygen.org/ln2data/releases/data/current/reviewed/"

GTEX_V8_MEDIAN_TPM_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
)

# GlyGen reviewed flat-file catalog: each value is (set_name, filename, description).
# Source: GlyGen, NIH/NIGMS Common Fund Glycoscience Program.
# https://data.glygen.org/ln2data/releases/data/current/reviewed/
GLYGEN_CATALOGS = [
    (
        "GlyGen_glycosyltransferases",
        "human_protein_glycosyltransferase.csv",
        (
            "Human glycosyltransferase genes from GlyGen reviewed data "
            "(NIH/NIGMS Common Fund Glycoscience Program; public). "
            "Enzymes that transfer sugar moieties to acceptor substrates. "
            "Source: https://data.glygen.org/ln2data/releases/data/current/reviewed/"
            "human_protein_glycosyltransferase.csv"
        ),
    ),
    (
        "GlyGen_glycohydrolases",
        "human_protein_glycohydrolase.csv",
        (
            "Human glycohydrolase (glycosidase) genes from GlyGen reviewed data "
            "(NIH/NIGMS Common Fund Glycoscience Program; public). "
            "Enzymes that hydrolyze glycosidic bonds. "
            "Source: https://data.glygen.org/ln2data/releases/data/current/reviewed/"
            "human_protein_glycohydrolase.csv"
        ),
    ),
    (
        "GlyGen_glycogenes",
        "human_protein_glycogenes.csv",
        (
            "Human glycogene catalogue from GlyGen reviewed data "
            "(NIH/NIGMS Common Fund Glycoscience Program; public). "
            "Comprehensive set of genes involved in glycan biosynthesis, modification, "
            "and recognition. "
            "Source: https://data.glygen.org/ln2data/releases/data/current/reviewed/"
            "human_protein_glycogenes.csv"
        ),
    ),
    (
        "GlyGen_glycosylation_motif_proteins",
        "human_protein_glycosylation_motifs.csv",
        (
            "Human proteins with defined glycosylation sequence motifs from GlyGen "
            "reviewed data (NIH/NIGMS Common Fund Glycoscience Program; public). "
            "Includes N-glycosylation sequons (Asn-X-Ser/Thr), O-glycosylation motifs, "
            "and GPI-anchor signal sequences. "
            "Source: https://data.glygen.org/ln2data/releases/data/current/reviewed/"
            "human_protein_glycosylation_motifs.csv"
        ),
    ),
]

# O-linked glycan catalogs derived from GlyGen annotation/xref files.
# Each entry: (set_name, description, fetch_fn_key)
# fetch_fn_key is used in run() to dispatch to the right specialized loader.
_O_LINKED_CATALOGS = [
    (
        "GlyGen_O_linked_glycoproteins",
        "o_linked_site_annotation",
        (
            "Human O-linked glycoprotein substrate genes from GlyGen reviewed data "
            "(NIH/NIGMS Common Fund Glycoscience Program; public). "
            "Proteins carrying experimentally confirmed O-linked glycans (GalNAc, GlcNAc, "
            "Fuc, Man, Xyl, Glc, or glycosaminoglycan chains) based on UniProtKB site "
            "annotations; excludes microbial-infection-mediated modifications. "
            "Source: https://data.glygen.org/ln2data/releases/data/current/reviewed/"
            "human_protein_site_annotation_uniprotkb.csv"
        ),
    ),
    (
        "GlyGen_O_GlcNAc_proteins",
        "oglcnac_atlas",
        (
            "Human O-GlcNAc-modified proteins from the O-GlcNAc Atlas "
            "(cross-referenced in GlyGen reviewed data; NIH/NIGMS Common Fund "
            "Glycoscience Program; public). "
            "Proteins with experimentally confirmed O-GlcNAc modification (O-linked "
            "N-acetylglucosamine on Ser/Thr residues). "
            "Source: https://data.glygen.org/ln2data/releases/data/current/reviewed/"
            "human_protein_xref_oglcnac_atlas.csv"
        ),
    ),
]

# Gene symbol column candidates, in priority order
_GENE_COL_CANDIDATES = ("gene_symbol", "gene_name", "gene", "hgnc_symbol")


def _detect_gene_col(fieldnames: list[str] | None) -> str | None:
    if not fieldnames:
        return None
    fl = {f.strip().lower(): f for f in fieldnames}
    for cand in _GENE_COL_CANDIDATES:
        if cand in fl:
            return fl[cand]
    return None


def _safe_name(s: str, max_len: int = 60) -> str:
    return "".join(c if c.isalnum() else "_" for c in s)[:max_len]


def _maybe_download(path_or_url: str, out_dir: Path) -> Path:
    if path_or_url.startswith(("http://", "https://")):
        dl_dir = out_dir / "downloads"
        dl_dir.mkdir(parents=True, exist_ok=True)
        dest = dl_dir / Path(path_or_url.rstrip("/").split("/")[-1])
        if not dest.exists():
            urllib.request.urlretrieve(path_or_url, dest)
        return dest
    return Path(path_or_url)


def _load_ac_to_gene(out_dir: Path) -> dict[str, str]:
    """Build uniprotkb_canonical_ac → gene_symbol_recommended map from GlyGen."""
    names_fname = "human_protein_genenames_uniprotkb.csv"
    local = _maybe_download(GLYGEN_BASE_URL + names_fname, out_dir)
    data = local.read_text(encoding="utf-8", errors="replace")
    ac_to_gene: dict[str, str] = {}
    for row in csv.DictReader(io.StringIO(data)):
        ac = (row.get("uniprotkb_canonical_ac") or "").strip()
        sym = (row.get("gene_symbol_recommended") or "").strip()
        if ac and sym and sym not in ("", "NA"):
            ac_to_gene.setdefault(ac, sym)
    return ac_to_gene


def _fetch_glygen_genes(fname: str, out_dir: Path) -> tuple[list[str], Path]:
    url = GLYGEN_BASE_URL + fname
    local = _maybe_download(url, out_dir)
    data = local.read_text(encoding="utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(data))
    gene_col = _detect_gene_col(reader.fieldnames)
    if gene_col is not None:
        genes = sorted({
            row[gene_col].strip()
            for row in reader
            if row.get(gene_col) and row[gene_col].strip() not in ("", "NA")
        })
    else:
        # No gene symbol column — join via uniprotkb_canonical_ac
        ac_to_gene = _load_ac_to_gene(out_dir)
        genes = sorted({
            ac_to_gene[row["uniprotkb_canonical_ac"].strip()]
            for row in reader
            if row.get("uniprotkb_canonical_ac")
            and row["uniprotkb_canonical_ac"].strip() in ac_to_gene
        })
    return genes, local


def _fetch_o_linked_glycoprotein_genes(out_dir: Path) -> tuple[list[str], Path]:
    """Genes carrying confirmed O-linked glycans from GlyGen site annotations.

    Filters human_protein_site_annotation_uniprotkb.csv for rows where
    annotation starts with "O-linked" (excludes microbial-infection prefix).
    """
    fname = "human_protein_site_annotation_uniprotkb.csv"
    url = GLYGEN_BASE_URL + fname
    local = _maybe_download(url, out_dir)
    data = local.read_text(encoding="utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(data))
    genes: set[str] = set()
    for row in reader:
        ann = (row.get("annotation") or "").strip()
        if not ann.startswith("O-linked"):
            continue
        gene = (row.get("gene_symbol") or "").strip()
        if gene and gene not in ("", "NA"):
            genes.add(gene)
    return sorted(genes), local


def _fetch_oglcnac_genes(out_dir: Path) -> tuple[list[str], Path]:
    """Genes with O-GlcNAc modification from the O-GlcNAc Atlas via GlyGen xref.

    Joins human_protein_xref_oglcnac_atlas.csv with
    human_protein_genenames_uniprotkb.csv on uniprotkb_canonical_ac to get
    gene symbols.
    """
    xref_fname = "human_protein_xref_oglcnac_atlas.csv"
    names_fname = "human_protein_genenames_uniprotkb.csv"
    xref_local = _maybe_download(GLYGEN_BASE_URL + xref_fname, out_dir)
    names_local = _maybe_download(GLYGEN_BASE_URL + names_fname, out_dir)

    # Build AC → gene symbol map from the gene names file
    ac_to_gene: dict[str, str] = {}
    names_data = names_local.read_text(encoding="utf-8", errors="replace")
    for row in csv.DictReader(io.StringIO(names_data)):
        ac = (row.get("uniprotkb_canonical_ac") or "").strip()
        sym = (row.get("gene_symbol_recommended") or "").strip()
        if ac and sym and sym not in ("", "NA"):
            ac_to_gene.setdefault(ac, sym)

    xref_data = xref_local.read_text(encoding="utf-8", errors="replace")
    genes: set[str] = set()
    for row in csv.DictReader(io.StringIO(xref_data)):
        ac = (row.get("uniprotkb_canonical_ac") or "").strip()
        sym = ac_to_gene.get(ac)
        if sym:
            genes.add(sym)
    return sorted(genes), xref_local


def _load_tstat_matrix(path: Path) -> tuple[list[str], list[str], dict[str, list[float]]]:
    """Load a gene × tissue numeric matrix (t-stats or median TPM).

    Accepts:
    - Plain TSV: gene<TAB>tissue1<TAB>tissue2...
    - GCT v1.2 (gzip or plain): Name<TAB>Description<TAB>tissue1... with two
      preamble lines (#1.2 and dimensions). The Description column is skipped
      and gene names are taken from the Name/first column.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header: list[str] = []
        for row in reader:
            if not row or not row[0]:
                continue
            # Skip GCT preamble: '#1.2' line and integer dimension line
            if row[0].startswith("#"):
                continue
            try:
                int(row[0])
                continue  # dimension line (e.g. "56200\t54")
            except ValueError:
                pass
            header = row
            break

        # Detect GCT format: second column is 'Description'/'description'.
        # In GTEx GCT files the Description column holds the gene symbol
        # (e.g. WASH7P), while the Name column holds the ENSG ID.
        # Use the Description value as the gene key so GlyGen symbol lookups work.
        gct = len(header) > 1 and header[1].lower() in ("description", "name", "id")
        val_start = 2 if gct else 1
        gene_col = 1 if gct else 0  # symbol is in Description for GCT
        tissues = header[val_start:]

        genes: list[str] = []
        gene_tstat: dict[str, list[float]] = {}
        for row in reader:
            if not row or not row[0]:
                continue
            gene = row[gene_col] if len(row) > gene_col else row[0]
            gene = gene.strip()
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


def _write_set(
    set_dir: Path,
    set_name: str,
    description: str,
    genes: list[str],
    converter_name: str,
    parameters: dict,
    files: list,
    set_subtype: str,
) -> None:
    set_dir.mkdir(parents=True, exist_ok=True)
    with (set_dir / "geneset.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\n")
        for g in genes:
            fh.write(f"{g}\n")
    write_gmt([(set_name, genes)], set_dir / "genesets.gmt")
    meta = make_metadata(
        converter_name=converter_name,
        parameters=parameters,
        data_type="glycoproteomics_catalog",
        assay="curated_database",
        organism="human",
        genome_build="GRCh38",
        files=files,
        gene_annotation={
            "mode": "curated",
            "source": "glygen_reviewed",
            "database_url": GLYGEN_BASE_URL,
            "funding": "NIH/NIGMS Common Fund Glycoscience Program",
        },
        weights={
            "weight_type": "nonnegative",
            "normalization": {"method": "none", "target_sum": None},
            "aggregation": "catalog_membership",
        },
        summary={
            "n_input_features": len(genes),
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
    activate_runtime_context("glygen_gtex", getattr(args, "provenance_overlay_json", None))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tstat_path = _maybe_download(args.gtex_tstat_tsv, out_dir)
    tissues, _all_genes, gene_tstat = _load_tstat_matrix(tstat_path)
    threshold = float(getattr(args, "tstat_threshold", 4.0))

    gtex_file_rec = input_file_record(str(tstat_path), "gtex_tstat_tsv")

    n_standalone = 0
    n_x_gtex = 0
    standalone_counts: list[int] = []
    x_gtex_counts: list[int] = []

    for set_name, fname, standalone_desc in GLYGEN_CATALOGS:
        genes, glygen_local = _fetch_glygen_genes(fname, out_dir)
        glygen_url = GLYGEN_BASE_URL + fname
        glygen_file_rec = input_file_record(str(glygen_local), "glygen_reviewed_csv")

        # Standalone set (no GTEx filter — all genes in this GlyGen category)
        standalone_dir = out_dir / "standalone" / set_name
        _write_set(
            set_dir=standalone_dir,
            set_name=set_name,
            description=standalone_desc,
            genes=genes,
            converter_name="glygen_gtex",
            parameters={
                "set_name": set_name,
                "glygen_file": fname,
                "glygen_base_url": GLYGEN_BASE_URL,
                "glygen_file_url": glygen_url,
                "set_subtype": "standalone",
                "funding": "NIH/NIGMS Common Fund Glycoscience Program (GlyGen)",
            },
            files=[glygen_file_rec],
            set_subtype="standalone",
        )
        standalone_counts.append(len(genes))
        n_standalone += 1

        # x GTEx: per-tissue subsets where t-stat >= threshold
        for ti, tissue in enumerate(tissues):
            enriched = sorted(
                g for g in genes
                if g in gene_tstat and gene_tstat[g][ti] >= threshold
            )
            if not enriched:
                continue
            x_name = f"{set_name}_x_GTEx_enriched_{_safe_name(tissue)}"
            x_desc = (
                f"{set_name} genes with GTEx tissue-enrichment (t-stat >= {threshold}) "
                f"in {tissue} (GTEx V8 median TPM; NIH Common Fund). "
                f"Derived from GlyGen reviewed file {fname} "
                f"(NIH/NIGMS Common Fund Glycoscience Program; "
                f"{glygen_url}) "
                f"intersected with GTEx V8 tissue specificity. "
                f"Absence from this set means t-stat < {threshold} in this tissue, "
                f"not that the gene is absent from the organism."
            )
            x_dir = out_dir / "x_gtex" / x_name
            _write_set(
                set_dir=x_dir,
                set_name=x_name,
                description=x_desc,
                genes=enriched,
                converter_name="glygen_gtex",
                parameters={
                    "set_name": set_name,
                    "glygen_file": fname,
                    "glygen_base_url": GLYGEN_BASE_URL,
                    "glygen_file_url": glygen_url,
                    "tissue": tissue,
                    "tstat_threshold": threshold,
                    "gtex_version": "v8",
                    "gtex_source_url": GTEX_V8_MEDIAN_TPM_URL,
                    "set_subtype": "x_gtex_tissue_enriched",
                    "funding": (
                        "NIH/NIGMS Common Fund Glycoscience Program (GlyGen) + "
                        "NIH Common Fund (GTEx)"
                    ),
                },
                files=[glygen_file_rec, gtex_file_rec],
                set_subtype="x_gtex_tissue_enriched",
            )
            x_gtex_counts.append(len(enriched))
            n_x_gtex += 1

    # O-linked glycan catalogs (specialized loaders)
    _o_linked_fetchers = {
        "o_linked_site_annotation": _fetch_o_linked_glycoprotein_genes,
        "oglcnac_atlas": _fetch_oglcnac_genes,
    }
    for set_name, fetch_key, standalone_desc in _O_LINKED_CATALOGS:
        fetch_fn = _o_linked_fetchers[fetch_key]
        genes, glygen_local = fetch_fn(out_dir)
        if not genes:
            continue
        # Derive source filename from the local file for provenance
        fname = glygen_local.name
        glygen_url = GLYGEN_BASE_URL + fname
        glygen_file_rec = input_file_record(str(glygen_local), "glygen_reviewed_csv")

        standalone_dir = out_dir / "standalone" / set_name
        _write_set(
            set_dir=standalone_dir,
            set_name=set_name,
            description=standalone_desc,
            genes=genes,
            converter_name="glygen_gtex",
            parameters={
                "set_name": set_name,
                "glygen_file": fname,
                "glygen_base_url": GLYGEN_BASE_URL,
                "glygen_file_url": glygen_url,
                "set_subtype": "standalone",
                "funding": "NIH/NIGMS Common Fund Glycoscience Program (GlyGen)",
            },
            files=[glygen_file_rec],
            set_subtype="standalone",
        )
        standalone_counts.append(len(genes))
        n_standalone += 1

        for ti, tissue in enumerate(tissues):
            enriched = sorted(
                g for g in genes
                if g in gene_tstat and gene_tstat[g][ti] >= threshold
            )
            if not enriched:
                continue
            x_name = f"{set_name}_x_GTEx_enriched_{_safe_name(tissue)}"
            x_desc = (
                f"{set_name} genes with GTEx tissue-enrichment (t-stat >= {threshold}) "
                f"in {tissue} (GTEx V8 median TPM; NIH Common Fund). "
                f"Derived from GlyGen reviewed file {fname} "
                f"(NIH/NIGMS Common Fund Glycoscience Program; "
                f"{glygen_url}) "
                f"intersected with GTEx V8 tissue specificity."
            )
            x_dir = out_dir / "x_gtex" / x_name
            _write_set(
                set_dir=x_dir,
                set_name=x_name,
                description=x_desc,
                genes=enriched,
                converter_name="glygen_gtex",
                parameters={
                    "set_name": set_name,
                    "glygen_file": fname,
                    "glygen_base_url": GLYGEN_BASE_URL,
                    "glygen_file_url": glygen_url,
                    "tissue": tissue,
                    "tstat_threshold": threshold,
                    "gtex_version": "v8",
                    "gtex_source_url": GTEX_V8_MEDIAN_TPM_URL,
                    "set_subtype": "x_gtex_tissue_enriched",
                    "funding": (
                        "NIH/NIGMS Common Fund Glycoscience Program (GlyGen) + "
                        "NIH Common Fund (GTEx)"
                    ),
                },
                files=[glygen_file_rec, gtex_file_rec],
                set_subtype="x_gtex_tissue_enriched",
            )
            x_gtex_counts.append(len(enriched))
            n_x_gtex += 1

    return {
        "n_standalone": n_standalone,
        "n_x_gtex": n_x_gtex,
        "n_groups": n_standalone + n_x_gtex,
        "n_genes_standalone_min": min(standalone_counts) if standalone_counts else 0,
        "n_genes_standalone_max": max(standalone_counts) if standalone_counts else 0,
        "n_genes_x_gtex_min": min(x_gtex_counts) if x_gtex_counts else 0,
        "n_genes_x_gtex_max": max(x_gtex_counts) if x_gtex_counts else 0,
        "out_dir": str(out_dir),
    }
