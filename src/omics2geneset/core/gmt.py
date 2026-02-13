from __future__ import annotations

from pathlib import Path
import re


_WS_RE = re.compile(r"\s+")
_BAD_RE = re.compile(r"[^A-Za-z0-9_.=:-]")
_ENSEMBL_LIKE_RE = re.compile(r"^ENS[A-Z]*G[0-9]+(?:\.[0-9]+)?$")


def sanitize_gmt_name(name: str) -> str:
    s = _WS_RE.sub("_", str(name).strip())
    s = s.replace("/", "_").replace("\\", "_")
    s = _BAD_RE.sub("", s)
    return s or "geneset"


def _is_ensembl_like(value: str) -> bool:
    return bool(_ENSEMBL_LIKE_RE.match(value))


def _candidate_token(row: dict[str, object], prefer_symbol: bool, require_symbol: bool = False) -> str:
    gene_id = str(row.get("gene_id", "")).strip()
    gene_symbol_obj = row.get("gene_symbol")
    gene_symbol = "" if gene_symbol_obj is None else str(gene_symbol_obj).strip()
    if require_symbol:
        if gene_symbol and not _is_ensembl_like(gene_symbol):
            return gene_symbol
        return ""
    if prefer_symbol and gene_symbol:
        return gene_symbol
    return gene_id


def choose_gene_tokens(rows, prefer_symbol: bool, require_symbol: bool = False) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for row in rows:
        token = _candidate_token(row, prefer_symbol, require_symbol=require_symbol)
        if not token or token in seen:
            continue
        seen.add(token)
        out.append(token)
    return out


def write_gmt(gene_sets: list[tuple[str, list[str]]], out_path: str | Path) -> None:
    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8", newline="") as fh:
        for name, genes in gene_sets:
            fh.write(f"{sanitize_gmt_name(name)}\t{' '.join(genes)}\n")


def parse_int_list_csv(value: str) -> list[int]:
    text = str(value).strip()
    if not text:
        return []
    out: list[int] = []
    for part in text.split(","):
        tok = part.strip()
        if not tok:
            continue
        parsed = int(tok)
        if parsed <= 0:
            raise ValueError("gmt_topk_list values must be positive integers")
        out.append(parsed)
    return out


def parse_mass_list_csv(value: str) -> list[float]:
    text = str(value).strip()
    if not text:
        return []
    out: list[float] = []
    for part in text.split(","):
        tok = part.strip()
        if not tok:
            continue
        parsed = float(tok)
        if not (0.0 < parsed < 1.0):
            raise ValueError("gmt_mass_list values must be in (0,1)")
        out.append(parsed)
    return out


def parse_str_list_csv(value: str) -> list[str]:
    text = str(value).strip()
    if not text:
        return []
    return [tok.strip() for tok in text.split(",") if tok.strip()]


def resolve_gmt_out_path(base_dir: Path, gmt_out: str | None) -> Path:
    if gmt_out is None or not str(gmt_out).strip():
        return base_dir / "genesets.gmt"
    p = Path(gmt_out)
    if p.is_absolute():
        return p
    return base_dir / p


def _clamp_k(k: int, min_genes: int, max_genes: int) -> int:
    return max(min_genes, min(max_genes, int(k)))


def _rows_sorted_by_score(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return sorted(rows, key=lambda r: (-float(r.get("score", 0.0)), str(r.get("gene_id", ""))))


def _dup_count(rows: list[dict[str, object]], prefer_symbol: bool, require_symbol: bool, emitted: list[str]) -> int:
    candidates = [tok for tok in (_candidate_token(r, prefer_symbol, require_symbol=require_symbol) for r in rows) if tok]
    return max(0, len(candidates) - len(emitted))


def _passes_biotype_filter(row: dict[str, object], allowed_biotypes: set[str] | None) -> bool:
    if not allowed_biotypes:
        return True
    b = row.get("gene_biotype")
    if b is None:
        # Backward compatibility for toy rows without biotype annotation.
        return True
    biotype = str(b).strip().lower()
    if not biotype:
        return True
    return biotype in allowed_biotypes


def _topk_plan_rows(rows: list[dict[str, object]], k: int) -> list[dict[str, object]]:
    if k <= 0:
        return []
    return rows[: min(k, len(rows))]


def _mass_plan_rows(
    rows: list[dict[str, object]],
    tau: float,
    min_genes: int,
    max_genes: int,
) -> tuple[list[dict[str, object]], int]:
    ranked = sorted(rows, key=lambda r: (-max(float(r.get("score", 0.0)), 0.0), str(r.get("gene_id", ""))))
    nonneg = [max(float(r.get("score", 0.0)), 0.0) for r in ranked]
    total = sum(nonneg)
    if total <= 0:
        k = min_genes
    else:
        cum = 0.0
        k = 0
        for val in nonneg:
            k += 1
            cum += val / total
            if cum >= tau:
                break
    k = _clamp_k(k, min_genes, max_genes)
    return ranked[: min(k, len(ranked))], k


def build_gmt_sets_from_rows(
    rows: list[dict[str, object]],
    base_name: str,
    prefer_symbol: bool,
    min_genes: int,
    max_genes: int,
    topk_list: list[int],
    mass_list: list[float],
    split_signed: bool,
    require_symbol: bool = False,
    allowed_biotypes: set[str] | None = None,
) -> tuple[list[tuple[str, list[str]]], list[dict[str, object]]]:
    if min_genes <= 0 or max_genes <= 0:
        raise ValueError("gmt_min_genes and gmt_max_genes must be positive")
    if min_genes > max_genes:
        raise ValueError("gmt_min_genes must be <= gmt_max_genes")

    ranked = _rows_sorted_by_score([r for r in rows if _passes_biotype_filter(r, allowed_biotypes)])
    has_negative = any(float(r.get("score", 0.0)) < 0.0 for r in ranked)

    variants: list[tuple[str, list[dict[str, object]]]]
    if split_signed and has_negative:
        pos_rows = _rows_sorted_by_score(
            [{**r, "score": max(float(r.get("score", 0.0)), 0.0)} for r in ranked]
        )
        neg_rows = _rows_sorted_by_score(
            [{**r, "score": max(-float(r.get("score", 0.0)), 0.0)} for r in ranked]
        )
        variants = [("__pos", pos_rows), ("__neg", neg_rows)]
    else:
        variants = [("", ranked)]

    sets: list[tuple[str, list[str]]] = []
    plans: list[dict[str, object]] = []
    seen_names: set[str] = set()

    for sign_suffix, variant_rows in variants:
        for requested_k in topk_list:
            k = _clamp_k(int(requested_k), min_genes, max_genes)
            selected_rows = _topk_plan_rows(variant_rows, k)
            genes = choose_gene_tokens(selected_rows, prefer_symbol, require_symbol=require_symbol)
            set_name = sanitize_gmt_name(f"{base_name}{sign_suffix}__topk={k}")
            if set_name in seen_names or not genes:
                continue
            seen_names.add(set_name)
            sets.append((set_name, genes))
            plans.append(
                {
                    "name": set_name,
                    "method": "topk",
                    "parameters": {"k": k},
                    "n_genes_emitted": len(genes),
                    "token_type": "gene_symbol" if prefer_symbol else "gene_id",
                    "n_duplicates_dropped": _dup_count(
                        selected_rows,
                        prefer_symbol,
                        require_symbol,
                        genes,
                    ),
                }
            )

        for tau in mass_list:
            selected_rows, k = _mass_plan_rows(variant_rows, float(tau), min_genes, max_genes)
            genes = choose_gene_tokens(selected_rows, prefer_symbol, require_symbol=require_symbol)
            tau_str = format(float(tau), ".6g")
            set_name = sanitize_gmt_name(f"{base_name}{sign_suffix}__hpd_mass={tau_str}__k={k}")
            if set_name in seen_names or not genes:
                continue
            seen_names.add(set_name)
            sets.append((set_name, genes))
            plans.append(
                {
                    "name": set_name,
                    "method": "hpd_mass",
                    "parameters": {"tau": float(tau), "k": k},
                    "n_genes_emitted": len(genes),
                    "token_type": "gene_symbol" if prefer_symbol else "gene_id",
                    "n_duplicates_dropped": _dup_count(
                        selected_rows,
                        prefer_symbol,
                        require_symbol,
                        genes,
                    ),
                }
            )
    return sets, plans
