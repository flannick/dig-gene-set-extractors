#!/usr/bin/env python3
"""KidsFirst curated disease gene-set workflow (DIG-owned).

Formalizes the KidsFirst / CBTN pediatric-cancer "curated disease-up" gene-set
selection into a single reusable DIG workflow: read the per-comparison
``deg_long.tsv`` differential-expression tables, take the concordant
upregulated set across a disease's controls, apply the mean-signed-neglog10padj
score threshold and safety cap, and emit the disease gene sets + auxiliary
tissue/subtype marker sets + a curation manifest. These curated products are the
KidsFirst library's HZ2 model outputs.

This logic was migrated verbatim from
geneset-extractor-dev/KidsFirst/src/curate_disease_genesets.py so that DIG owns the
reusable KidsFirst workflow and geneset-extractor-dev/KidsFirst acts as the
config/wrapper layer (branch two-repo standard) — the same migration already done
for ``kidsfirst_prepare``.

Invoke via the DIG CLI:
    geneset-extractors workflows kidsfirst_curate \
        --analysis_dir <dir with {COMPARISON}/de_results/deg_long.tsv> \
        --out_dir <out> [--score_threshold 2.0 --safety_cap 200 ...]
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

from geneset_extractors.workflows.gtex_runtime_common import write_workflow_provenance_graph


# ── Thresholds ───────────────────────────────────────────────────────────────
PADJ_MAX_DEFAULT = 0.05
MIN_LOGFC_DEFAULT = 1.0
SCORE_THRESHOLD_DEFAULT = 2.0    # mean signed_neglog10padj >= this (≈ padj < 0.01)
                                  # applied after concordance; both filters together
SCORE_THRESHOLD_RELAXED = 1.30103  # fallback if result < MIN_GENES (≈ padj < 0.05)
SAFETY_CAP_DEFAULT = 200         # ceiling only — not a target
MIN_GENES_DEFAULT = 50           # warn + relax threshold if set is smaller than this
MIN_CORE_DEFAULT = 10            # fall back from intersection if concordant set < this


# ── Disease config ────────────────────────────────────────────────────────────
# compare_type:
#   "tumor_vs_normal"   → down genes = tissue/control markers, NOT disease biology
#   "subtype_contrast"  → down genes may carry biologically meaningful signal
#
# Multiple comparisons → concordant up-set (intersection/majority/union).
# Missing DE results are skipped gracefully; run sbatch_04 again after sbatch_02.
DISEASE_CONFIG = [
    # ── Kids First: blood / lymphoid ─────────────────────────────────────────
    {
        "label": "KF_TALL",
        "compare_type": "tumor_vs_normal",
        "comparisons": ["KF-TALL-vs-T21", "KF-TALL-vs-GTEx"],
        "note": "T-ALL concordant upregulated (T21 pediatric + GTEx whole-blood controls)",
    },
    # ── Kids First: solid tumors ──────────────────────────────────────────────
    {
        "label": "KF_NBL",
        "compare_type": "tumor_vs_normal",
        "comparisons": ["KF-NBL-vs-adrenal"],
        "note": "Neuroblastoma vs GTEx adrenal gland",
    },
    {
        "label": "KF_ESGR",
        "compare_type": "tumor_vs_normal",
        "comparisons": ["KF-ESGR-vs-muscle"],
        "note": "Ewing sarcoma vs GTEx skeletal muscle",
    },
    {
        "label": "KF_MMC",
        "compare_type": "tumor_vs_normal",
        "comparisons": ["KF-MMC-vs-blood"],
        "note": "AML vs GTEx whole blood",
    },
    # ── CBTN: brain tumors (each vs GTEx brain cortex) ────────────────────────
    # Auto-skipped until sbatch_02_cbtn_de.sh completes.
    # Results land in the same ANALYSIS_DIR as KidsFirst comparisons.
    {
        "label": "CBTN_LGG",
        "compare_type": "tumor_vs_normal",
        "comparisons": ["CBTN-low_grade_glioma-vs-brain_cortex"],
        "note": "Low-grade glioma vs GTEx brain cortex",
    },
    {
        "label": "CBTN_MG",
        "compare_type": "tumor_vs_normal",
        "comparisons": ["CBTN-malignant_glioma-vs-brain_cortex"],
        "note": "Malignant glioma vs GTEx brain cortex",
    },
    {
        "label": "CBTN_MB",
        "compare_type": "tumor_vs_normal",
        "comparisons": ["CBTN-medulloblastoma-vs-brain_cortex"],
        "note": "Medulloblastoma vs GTEx brain cortex",
    },
    {
        "label": "CBTN_EP",
        "compare_type": "tumor_vs_normal",
        "comparisons": ["CBTN-ependymoma-vs-brain_cortex"],
        "note": "Ependymoma vs GTEx brain cortex",
    },
    {
        "label": "CBTN_GG",
        "compare_type": "tumor_vs_normal",
        "comparisons": ["CBTN-ganglioglioma-vs-brain_cortex"],
        "note": "Ganglioglioma vs GTEx brain cortex",
    },
    {
        "label": "CBTN_CP",
        "compare_type": "tumor_vs_normal",
        "comparisons": ["CBTN-craniopharyngioma-vs-brain_cortex"],
        "note": "Craniopharyngioma vs GTEx brain cortex",
    },
    {
        "label": "CBTN_ATRT",
        "compare_type": "tumor_vs_normal",
        "comparisons": ["CBTN-atypical_teratoid_rhabdoid_tumor-vs-brain_cortex"],
        "note": "ATRT vs GTEx brain cortex",
    },
]


# ── Helpers ───────────────────────────────────────────────────────────────────

def _signed_neglog10padj(logfc: float, padj: float, floor: float = 1e-300) -> float:
    return math.copysign(-math.log10(max(padj, floor)), logfc)


def _is_valid_symbol(sym: str) -> bool:
    return bool(sym) and sym.upper() not in ("NA", "NULL", "NONE", "")


def _load_deg(
    path: Path,
    padj_max: float,
    min_logfc: float,
) -> tuple[dict[str, dict], dict[str, dict]]:
    """
    Parse deg_long.tsv (columns: comparison_id, gene_id, gene_symbol, logFC,
    stat, pvalue, padj). Returns (up_dict, dn_dict) keyed by gene_symbol.
    For duplicate symbols, keeps the entry with highest |score|.
    """
    up: dict[str, dict] = {}
    dn: dict[str, dict] = {}

    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sym = row.get("gene_symbol", "").strip()
            if not _is_valid_symbol(sym):
                continue
            try:
                logfc = float(row["logFC"])
                padj_raw = row.get("padj", "").strip()
                if not padj_raw or padj_raw == "NA":
                    continue
                padj = float(padj_raw)
            except (ValueError, KeyError):
                continue

            if padj >= padj_max:
                continue

            score = _signed_neglog10padj(logfc, padj)
            entry = {"logfc": logfc, "padj": padj, "score": score}

            if logfc >= min_logfc:
                existing = up.get(sym)
                if existing is None or abs(score) > abs(existing["score"]):
                    up[sym] = entry
            elif logfc <= -min_logfc:
                existing = dn.get(sym)
                if existing is None or abs(score) > abs(existing["score"]):
                    dn[sym] = entry

    return up, dn


def _concordant_set(
    per_comp: list[dict[str, dict]],
    min_core: int,
) -> tuple[set[str], str]:
    """
    Concordant gene set across comparisons.
    Strategy order: intersection → majority vote → union.
    """
    n = len(per_comp)
    if n == 1:
        return set(per_comp[0].keys()), "single"

    sets = [set(d.keys()) for d in per_comp]

    core = sets[0].intersection(*sets[1:])
    if len(core) >= min_core:
        return core, "intersection"

    votes: dict[str, int] = defaultdict(int)
    for s in sets:
        for sym in s:
            votes[sym] += 1
    thresh = n // 2 + 1
    majority = {sym for sym, cnt in votes.items() if cnt >= thresh}
    if len(majority) >= min_core:
        return majority, f"majority(>={thresh}/{n})"

    return set().union(*sets), "union"


def _rank_genes(
    gene_set: set[str],
    per_comp: list[dict[str, dict]],
    score_threshold: float,
    safety_cap: int,
    direction: str = "up",
) -> tuple[list[tuple[str, float]], int]:
    """
    Score each gene as mean signed_neglog10padj across controls where it appears.
    up → sort descending; down → sort ascending (most negative first).
    Apply score_threshold filter, then safety_cap.
    Returns (ranked_list, n_before_cap).
    """
    scored = []
    for sym in gene_set:
        scores = [d[sym]["score"] for d in per_comp if sym in d]
        if not scores:
            continue
        mean_score = sum(scores) / len(scores)
        if direction == "up" and mean_score >= score_threshold:
            scored.append((sym, mean_score))
        elif direction == "down" and mean_score <= -score_threshold:
            scored.append((sym, mean_score))

    scored.sort(key=lambda x: x[1], reverse=(direction == "up"))
    n_before_cap = len(scored)
    return scored[:safety_cap], n_before_cap


def _write_gmt_line(fh, set_id: str, description: str, genes: list[str]) -> None:
    fh.write("\t".join([set_id, description] + genes) + "\n")


def _write_tsv(path: Path, genes: list[tuple[str, float]], set_id: str, description: str) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene_symbol", "mean_score", "gene_set_id", "description"])
        for sym, score in genes:
            writer.writerow([sym, f"{score:.4f}", set_id, description])


# ── Per-disease curation ──────────────────────────────────────────────────────

def curate_disease(
    disease: dict,
    analysis_dir: Path,
    out_dir: Path,
    padj_max: float,
    min_logfc: float,
    score_threshold: float,
    safety_cap: int,
    min_genes: int,
    min_core: int,
) -> dict:
    label = disease["label"]
    compare_type = disease.get("compare_type", "tumor_vs_normal")
    comparisons = disease["comparisons"]
    note = disease.get("note", "")

    print(f"\n[{label}]  ({compare_type})", file=sys.stderr)

    up_per_comp: list[dict] = []
    dn_per_comp: list[dict] = []
    missing: list[str] = []

    for comp in comparisons:
        deg_path = analysis_dir / comp / "de_results" / "deg_long.tsv"
        if not deg_path.exists():
            print(f"  SKIP (missing): {deg_path}", file=sys.stderr)
            missing.append(comp)
            continue
        up, dn = _load_deg(deg_path, padj_max, min_logfc)
        print(f"  {comp}: {len(up)} up / {len(dn)} down (valid symbols)", file=sys.stderr)
        up_per_comp.append(up)
        dn_per_comp.append(dn)

    if not up_per_comp:
        return {"label": label, "compare_type": compare_type, "status": "missing", "missing": missing}

    # Disease gene set: concordance + score threshold
    up_set, up_strategy = _concordant_set(up_per_comp, min_core)

    up_ranked, n_before_cap = _rank_genes(up_set, up_per_comp, score_threshold, safety_cap, direction="up")
    up_threshold_used = score_threshold

    # Fallback: relax score threshold if result is too small
    if len(up_ranked) < min_genes:
        print(
            f"  [FALLBACK] {label}: {len(up_ranked)} genes at score>={score_threshold:.2f}"
            f" — relaxing to score>={SCORE_THRESHOLD_RELAXED:.5f} (padj<{padj_max})",
            file=sys.stderr,
        )
        up_ranked, n_before_cap = _rank_genes(
            up_set, up_per_comp, SCORE_THRESHOLD_RELAXED, safety_cap, direction="up"
        )
        up_threshold_used = SCORE_THRESHOLD_RELAXED

    n_up = len(up_ranked)
    warn_small = n_up < min_genes
    if warn_small:
        print(
            f"  [WARN] {label}: only {n_up} genes after fallback"
            f" — check gene symbol coverage",
            file=sys.stderr,
        )

    # Auxiliary down set (same logic, QC only for tumor_vs_normal)
    dn_set, dn_strategy = _concordant_set(dn_per_comp, min_core)
    dn_ranked, _ = _rank_genes(dn_set, dn_per_comp, score_threshold, safety_cap, direction="down")
    if len(dn_ranked) < min_genes:
        dn_ranked, _ = _rank_genes(
            dn_set, dn_per_comp, SCORE_THRESHOLD_RELAXED, safety_cap, direction="down"
        )
    n_dn = len(dn_ranked)

    # Write per-disease TSVs
    up_id = f"{label}_disease_up"
    dn_id = f"{label}_{'tissue_dn' if compare_type == 'tumor_vs_normal' else 'subtype_dn'}"
    _write_tsv(out_dir / f"{up_id}.tsv", up_ranked, up_id, note)
    _write_tsv(out_dir / f"{dn_id}.tsv", dn_ranked, dn_id, note)

    capped = n_before_cap > safety_cap
    print(
        f"  up  [{up_strategy}]: {len(up_set)} concordant"
        f" → score>={up_threshold_used:.2f}: {n_before_cap}"
        f" → final: {n_up}"
        + (" (capped)" if capped else "")
        + (" *** <50 ***" if warn_small else ""),
        file=sys.stderr,
    )
    print(
        f"  dn  [{dn_strategy}]: {len(dn_set)} concordant → {n_dn} final"
        + ("  (QC: tissue markers)" if compare_type == "tumor_vs_normal" else ""),
        file=sys.stderr,
    )

    return {
        "label": label,
        "compare_type": compare_type,
        "status": "ok",
        "n_comparisons_found": len(up_per_comp),
        "missing": missing,
        "warn_small": warn_small,
        "up_strategy": up_strategy,
        "up_threshold_used": up_threshold_used,
        "n_up_concordant": len(up_set),
        "n_up_score_pass": n_before_cap,
        "n_up_final": n_up,
        "up_genes": [sym for sym, _ in up_ranked],
        "dn_strategy": dn_strategy,
        "n_dn_concordant": len(dn_set),
        "n_dn_final": n_dn,
        "dn_genes": [sym for sym, _ in dn_ranked],
        "note": note,
    }


# ── Workflow entry point ───────────────────────────────────────────────────────

def run(args: argparse.Namespace) -> dict:
    analysis_dir = Path(args.analysis_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for disease in DISEASE_CONFIG:
        r = curate_disease(
            disease=disease,
            analysis_dir=analysis_dir,
            out_dir=out_dir,
            padj_max=args.padj_max,
            min_logfc=args.min_logfc,
            score_threshold=args.score_threshold,
            safety_cap=args.safety_cap,
            min_genes=args.min_genes,
            min_core=args.min_core,
        )
        results.append(r)

    # Write GMT files
    up_gmt    = out_dir / "disease_up.gmt"           # PRIMARY: disease gene sets
    tvn_dn    = out_dir / "tissue_markers_dn.gmt"    # QC: tumor_vs_normal down genes
    sub_dn    = out_dir / "subtype_markers_dn.gmt"   # auxiliary: subtype down genes

    with open(up_gmt, "w") as up_fh, \
         open(tvn_dn, "w") as tvn_fh, \
         open(sub_dn, "w") as sub_fh:
        for r in results:
            if r["status"] != "ok":
                continue
            label = r["label"]
            note = r.get("note", "")
            ctype = r["compare_type"]

            up_desc = (
                f"{note} | {r['up_strategy']}"
                f" concordant={r['n_up_concordant']} final={r['n_up_final']}"
            )
            dn_label = "tissue_dn" if ctype == "tumor_vs_normal" else "subtype_dn"
            dn_type_note = "(QC-tissue-markers)" if ctype == "tumor_vs_normal" else "(subtype-contrast)"
            dn_desc = (
                f"{note} {dn_type_note} | {r['dn_strategy']}"
                f" concordant={r['n_dn_concordant']} final={r['n_dn_final']}"
            )

            if r["up_genes"]:
                _write_gmt_line(up_fh, f"{label}_disease_up", up_desc, r["up_genes"])

            dn_target = tvn_fh if ctype == "tumor_vs_normal" else sub_fh
            if r["dn_genes"]:
                _write_gmt_line(dn_target, f"{label}_{dn_label}", dn_desc, r["dn_genes"])

    # Write manifest
    fieldnames = [
        "label", "compare_type", "status", "n_comparisons_found", "missing",
        "warn_small", "up_threshold_used",
        "up_strategy", "n_up_concordant", "n_up_score_pass", "n_up_final",
        "dn_strategy", "n_dn_concordant", "n_dn_final",
        "note",
    ]
    with open(out_dir / "manifest.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames,
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for r in results:
            row = {k: r.get(k, "") for k in fieldnames}
            if isinstance(row.get("missing"), list):
                row["missing"] = ",".join(row["missing"])
            writer.writerow(row)

    # Summary
    ok     = [r for r in results if r["status"] == "ok"]
    skip   = [r for r in results if r["status"] == "missing"]
    warn   = [r for r in ok if r.get("warn_small")]

    print("\n====== Curation summary ======", file=sys.stderr)
    print(
        f"  padj_max={args.padj_max}  min_logfc={args.min_logfc}"
        f"  score_threshold={args.score_threshold}  safety_cap={args.safety_cap}"
        f"  min_genes={args.min_genes}",
        file=sys.stderr,
    )
    print(
        f"\n  {'Disease':<20} {'Strategy':<20} {'Concord':>8} {'ScorePass':>10} {'Final':>6} {'Flag'}",
        file=sys.stderr,
    )
    print("  " + "-" * 78, file=sys.stderr)
    for r in results:
        if r["status"] == "missing":
            print(f"  {r['label']:<20} {'MISSING (skipped)'}", file=sys.stderr)
        else:
            flag = "*** <50 ***" if r.get("warn_small") else ""
            thr = r.get("up_threshold_used", "?")
            relaxed = " (relaxed)" if thr == SCORE_THRESHOLD_RELAXED else ""
            print(
                f"  {r['label']:<20} {r.get('up_strategy','?'):<20}"
                f" {r.get('n_up_concordant',0):>8}"
                f" {r.get('n_up_score_pass',0):>10}"
                f" {r.get('n_up_final',0):>6}"
                f"  {flag}{relaxed}",
                file=sys.stderr,
            )

    print(f"\n  {len(ok)} processed, {len(skip)} skipped (CBTN pending?), {len(warn)} with <50 up-genes", file=sys.stderr)
    print(f"\n  Outputs: {out_dir}", file=sys.stderr)
    print(f"    disease_up.gmt         PRIMARY delivery", file=sys.stderr)
    print(f"    tissue_markers_dn.gmt  QC only (tumor vs normal down = tissue markers)", file=sys.stderr)
    print(f"    subtype_markers_dn.gmt auxiliary (subtype contrast down, if any)", file=sys.stderr)
    print(f"    manifest.tsv           audit trail", file=sys.stderr)

    # Emit a native provenance graph for the curation step — real command, git-SHA
    # version, and md5/size on every deg_long input + curated output — so the shipped
    # HZ2 provenance is born clean (same mechanism as kidsfirst_prepare / the accepted
    # reference), instead of being hand-synthesized after the fact.
    deg_inputs: list[tuple[Path, str]] = []
    seen_paths: set[str] = set()
    for disease in DISEASE_CONFIG:
        for comp in disease["comparisons"]:
            p = (analysis_dir / comp / "de_results" / "deg_long.tsv").resolve()
            if p.exists() and str(p) not in seen_paths:
                seen_paths.add(str(p))
                deg_inputs.append((p, f"deg_long:{comp}"))
    curated_outputs = [
        (up_gmt, "disease_up_gmt"),
        (tvn_dn, "tissue_markers_dn_gmt"),
        (sub_dn, "subtype_markers_dn_gmt"),
        (out_dir / "manifest.tsv", "curation_manifest"),
    ]
    curated_outputs = [(p, role) for p, role in curated_outputs if p.exists()]
    write_workflow_provenance_graph(
        workflow_name="kidsfirst_curate",
        module_name="geneset_extractors.workflows.kidsfirst_curate",
        output_dir=out_dir,
        focus_output_path=up_gmt,
        output_paths=curated_outputs,
        input_paths=deg_inputs,
        parameters={
            "score_threshold": args.score_threshold,
            "safety_cap": args.safety_cap,
            "padj_max": args.padj_max,
            "min_logfc": args.min_logfc,
            "min_genes": args.min_genes,
            "concordance": "per_disease",
        },
        analysis_description=(
            "Analysis step that curates concordant disease-up gene sets from per-comparison "
            "DE results (kidsfirst_curate) and emits disease_up.gmt."
        ),
    )
    return {
        "out_dir": str(out_dir),
        "n_ok": len(ok),
        "n_skipped": len(skip),
        "n_warn_small": len(warn),
        "diseases": [r["label"] for r in ok],
        "provenance_graph": str(out_dir / "disease_up.provenance_graph.json"),
    }


def add_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--analysis_dir", required=True, type=Path,
                        help="Directory with {COMPARISON}/de_results/deg_long.tsv (KF + CBTN share this dir)")
    parser.add_argument("--out_dir",         required=True, type=Path)
    parser.add_argument("--padj_max",        type=float, default=PADJ_MAX_DEFAULT)
    parser.add_argument("--min_logfc",       type=float, default=MIN_LOGFC_DEFAULT)
    parser.add_argument("--score_threshold", type=float, default=SCORE_THRESHOLD_DEFAULT,
                        help="Min mean signed_neglog10padj after concordance (default 2.0 ≈ padj<0.01). "
                             "Relaxed to SCORE_THRESHOLD_RELAXED if result < --min_genes.")
    parser.add_argument("--safety_cap",      type=int,   default=SAFETY_CAP_DEFAULT,
                        help="Ceiling on final gene set size (not a target).")
    parser.add_argument("--min_genes",       type=int,   default=MIN_GENES_DEFAULT,
                        help="Warn + relax score threshold if final set is smaller than this.")
    parser.add_argument("--min_core",        type=int,   default=MIN_CORE_DEFAULT,
                        help="Min genes in intersection before falling back to majority/union")


def main() -> int:
    parser = argparse.ArgumentParser(description="Curate disease gene sets from DE results.")
    add_flags(parser)
    result = run(parser.parse_args())
    print(
        f"kidsfirst_curate_completed out={result['out_dir']} "
        f"n_ok={result['n_ok']} n_skipped={result['n_skipped']} n_warn_small={result['n_warn_small']}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
