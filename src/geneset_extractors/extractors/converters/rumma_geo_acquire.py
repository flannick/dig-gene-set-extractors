"""Acquire pinned RummaGEO GraphQL records for offline selection."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from geneset_extractors.extractors.converters.rumma_geo_selection import GENE_TERMS

ENDPOINT = "https://rummageo.com/graphql"
NOTEBOOK_COMMIT = "965d3a7299cdeaa8d54740b31093b80cebd5523b"
QUERY = """query TermSearch($terms:[String]!,$offset:Int=0,$first:Int=10000){geneSetTermSearch(terms:$terms,offset:$offset,first:$first){nodes{id term gse platform pmid publishedDate sampleGroups title geneSetById{nGeneIds species}} totalCount}}"""


def _request(payload: dict[str, object], endpoint: str) -> dict[str, object]:
    request = Request(endpoint, data=json.dumps(payload).encode(), headers={"Content-Type": "application/json"})
    try:
        with urlopen(request, timeout=120) as response:
            result = json.loads(response.read().decode())
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"RummaGEO GraphQL request failed: {exc}") from exc
    if result.get("errors"):
        raise RuntimeError(f"RummaGEO GraphQL returned errors: {result['errors']}")
    return result


def acquire(terms: list[str], *, first: int, endpoint: str = ENDPOINT, request_fn=_request) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    records: list[dict[str, object]] = []; pages: list[dict[str, object]] = []
    for term in terms:
        offset = 0; total = None; retrieved = 0; requests = 0
        while total is None or offset < total:
            payload = request_fn({"query": QUERY, "variables": {"terms": [term], "offset": offset, "first": first}}, endpoint)
            if payload.get("errors"):
                raise RuntimeError(f"RummaGEO GraphQL returned errors for {term!r}: {payload['errors']}")
            try: result = payload["data"]["geneSetTermSearch"]; nodes = result["nodes"]; total = int(result["totalCount"])
            except (KeyError, TypeError, ValueError) as exc: raise RuntimeError(f"invalid RummaGEO GraphQL response for {term!r}") from exc
            if not isinstance(nodes, list): raise RuntimeError(f"invalid nodes response for {term!r}")
            requests += 1
            for node in nodes: records.append({**node, "search_term": term})
            retrieved += len(nodes); offset += len(nodes)
            if not nodes and offset < total: raise RuntimeError(f"retrieved {retrieved} of {total} records for {term!r}")
        if retrieved != total: raise RuntimeError(f"retrieved {retrieved} records but API reported {total} for {term!r}")
        pages.append({"search_term": term, "total_count": total, "records_retrieved": retrieved, "multiple_pages": requests > 1})
    return sorted(records, key=lambda row: (str(row["search_term"]), str(row.get("id", "")), str(row.get("term", "")))), pages


def run(args):
    if args.model_id == "HZ2": terms = list(GENE_TERMS); drug_sha = None
    else:
        if not args.drug_terms_json: raise ValueError("HZ1 requires --drug_terms_json")
        drug_path = Path(args.drug_terms_json).resolve(); terms = json.loads(drug_path.read_text(encoding="utf-8")); drug_sha = hashlib.sha256(drug_path.read_bytes()).hexdigest()
    if not all(isinstance(term, str) for term in terms): raise ValueError("query terms must be strings")
    records, pages = acquire(terms, first=args.page_size, endpoint=args.endpoint)
    out = Path(args.out_dir).resolve(); out.mkdir(parents=True, exist_ok=True)
    cache = out / "rummageo_query_records.json"; cache.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    metadata = {"endpoint": args.endpoint, "retrieved_at": datetime.now(timezone.utc).isoformat(), "notebook_repository": "HarmonizomePythonScripts/RummaGEO", "notebook_commit": NOTEBOOK_COMMIT, "model_id": args.model_id, "search_terms": terms if args.model_id == "HZ2" else None, "drug_terms_sha256": drug_sha, "n_search_terms": len(terms), "n_query_records": len(records), "n_unique_uuids": len({str(row.get("id", "")) for row in records}), "pages": pages, "query_records_sha256": hashlib.sha256(cache.read_bytes()).hexdigest()}
    (out / "rummageo_query_records.acquisition.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"n_peaks": len(records), "n_genes": 0, "out_dir": str(out)}
