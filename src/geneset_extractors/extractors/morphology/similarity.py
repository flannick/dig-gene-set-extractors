from __future__ import annotations

import math


def align_feature_space(
    query_profiles: dict[str, dict[str, float]],
    reference_profiles: dict[str, dict[str, float]],
    *,
    feature_schema: list[str] | None = None,
) -> tuple[list[str], dict[str, object]]:
    query_features = set()
    for feat_map in query_profiles.values():
        query_features.update(feat_map)
    ref_features = set(feature_schema or [])
    if not ref_features:
        for feat_map in reference_profiles.values():
            ref_features.update(feat_map)
    shared = sorted(query_features & ref_features)
    frac = float(len(shared)) / float(len(ref_features)) if ref_features else 1.0
    return shared, {
        "n_query_features": len(query_features),
        "n_reference_features": len(ref_features),
        "n_shared_features": len(shared),
        "fraction_reference_features_shared": frac,
    }


def standardize_profiles(
    profiles: dict[str, dict[str, float]],
    *,
    features: list[str],
    feature_stats: dict[str, tuple[float, float]] | None,
) -> tuple[dict[str, list[float]], dict[str, object]]:
    out: dict[str, list[float]] = {}
    applied = feature_stats is not None
    for sample_id, feat_map in profiles.items():
        vec: list[float] = []
        for feat in features:
            val = float(feat_map.get(feat, 0.0))
            if feature_stats is not None and feat in feature_stats:
                center, scale = feature_stats[feat]
                denom = scale if scale != 0.0 else 1.0
                val = (val - center) / denom
            vec.append(float(val))
        out[sample_id] = vec
    return out, {"standardized": applied, "n_features": len(features)}


def cosine_similarity(a: list[float], b: list[float]) -> float:
    dot = sum(float(x) * float(y) for x, y in zip(a, b))
    na = math.sqrt(sum(float(x) * float(x) for x in a))
    nb = math.sqrt(sum(float(y) * float(y) for y in b))
    if na <= 0.0 or nb <= 0.0:
        return 0.0
    return float(dot / (na * nb))


def pearson_similarity(a: list[float], b: list[float]) -> float:
    if not a or not b or len(a) != len(b):
        return 0.0
    ma = sum(a) / float(len(a))
    mb = sum(b) / float(len(b))
    ca = [float(x) - ma for x in a]
    cb = [float(y) - mb for y in b]
    return cosine_similarity(ca, cb)


def pairwise_similarity(
    query_vectors: dict[str, list[float]],
    reference_vectors: dict[str, list[float]],
    *,
    metric: str,
) -> dict[str, dict[str, float]]:
    sim_fn = cosine_similarity if metric == "cosine" else pearson_similarity
    out: dict[str, dict[str, float]] = {}
    for query_id, qvec in query_vectors.items():
        out[query_id] = {}
        for ref_id, rvec in reference_vectors.items():
            out[query_id][ref_id] = float(sim_fn(qvec, rvec))
    return out
