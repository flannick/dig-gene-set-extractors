from __future__ import annotations

from geneset_extractors.extractors.converters.rumma_geo_acquire import acquire


def test_acquisition_paginates_and_preserves_search_term_uuid_relationship() -> None:
    calls = []
    def fake(payload, endpoint):
        calls.append(payload["variables"])
        term, offset = payload["variables"]["terms"][0], payload["variables"]["offset"]
        nodes = [{"id": "same", "term": f"{term}-{offset}"}] if offset < 2 else []
        return {"data": {"geneSetTermSearch": {"nodes": nodes, "totalCount": 2}}}
    records, pages = acquire(["b", "a"], first=1, request_fn=fake)
    assert [(row["search_term"], row["id"]) for row in records] == [("a", "same"), ("a", "same"), ("b", "same"), ("b", "same")]
    assert all(page["multiple_pages"] for page in pages)


def test_acquisition_fails_for_graphql_errors() -> None:
    def fake(payload, endpoint): return {"errors": [{"message": "bad query"}]}
    try: acquire(["knockdown"], first=10, request_fn=fake)
    except RuntimeError as exc: assert "errors" in str(exc)
    else: raise AssertionError("expected GraphQL failure")


def test_acquisition_detects_reported_count_mismatch() -> None:
    def fake(payload, endpoint): return {"data": {"geneSetTermSearch": {"nodes": [], "totalCount": 1}}}
    try: acquire(["knockdown"], first=10, request_fn=fake)
    except RuntimeError as exc: assert "retrieved 0 of 1" in str(exc)
    else: raise AssertionError("expected count mismatch")
