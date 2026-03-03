from __future__ import annotations


PROGRAM_PROMOTER_ACTIVITY = "promoter_activity"
PROGRAM_DISTAL_ACTIVITY = "distal_activity"
PROGRAM_LINKED_ACTIVITY = "linked_activity"

PROGRAM_PRESET_NONE = "none"
PROGRAM_PRESET_DEFAULT = "default"
PROGRAM_PRESET_CONNECTABLE = "connectable"
PROGRAM_PRESET_QC = "qc"
PROGRAM_PRESET_EXPERIMENTAL = "experimental"
PROGRAM_PRESET_ALL = "all"

ALLOWED_PROGRAM_METHODS = (
    PROGRAM_PROMOTER_ACTIVITY,
    PROGRAM_DISTAL_ACTIVITY,
    PROGRAM_LINKED_ACTIVITY,
)

_LINK_METHODS = (
    "promoter_overlap",
    "nearest_tss",
    "distance_decay",
)
_EXTERNAL_LINK_METHOD = "external"

_PRESETS: dict[str, tuple[str, ...]] = {
    PROGRAM_PRESET_NONE: (PROGRAM_LINKED_ACTIVITY,),
    PROGRAM_PRESET_DEFAULT: (
        PROGRAM_PROMOTER_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
    ),
    PROGRAM_PRESET_CONNECTABLE: (
        PROGRAM_PROMOTER_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
    ),
    PROGRAM_PRESET_QC: (
        PROGRAM_PROMOTER_ACTIVITY,
        PROGRAM_LINKED_ACTIVITY,
    ),
    PROGRAM_PRESET_EXPERIMENTAL: (
        PROGRAM_PROMOTER_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
        PROGRAM_LINKED_ACTIVITY,
    ),
    PROGRAM_PRESET_ALL: (
        PROGRAM_PROMOTER_ACTIVITY,
        PROGRAM_DISTAL_ACTIVITY,
        PROGRAM_LINKED_ACTIVITY,
    ),
}


def normalize_program_preset(preset: str | None) -> str:
    name = PROGRAM_PRESET_DEFAULT if preset is None else str(preset).strip()
    if not name:
        return PROGRAM_PRESET_CONNECTABLE
    if name == PROGRAM_PRESET_DEFAULT:
        return PROGRAM_PRESET_CONNECTABLE
    if name not in _PRESETS:
        raise ValueError(f"Unsupported program_preset: {name}")
    return name


def parse_csv_tokens(value: str | None) -> list[str]:
    text = "" if value is None else str(value).strip()
    if not text:
        return []
    return [tok.strip() for tok in text.split(",") if tok.strip()]


def resolve_program_methods(
    program_preset: str | None,
    program_methods_csv: str | None,
) -> list[str]:
    preset_name = normalize_program_preset(program_preset)
    requested = parse_csv_tokens(program_methods_csv)

    expanded: list[str] = []
    if requested:
        for token in requested:
            if token == "all":
                expanded.extend(ALLOWED_PROGRAM_METHODS)
            else:
                expanded.append(token)
    else:
        expanded.extend(_PRESETS[preset_name])

    out: list[str] = []
    seen: set[str] = set()
    for method in expanded:
        if method not in ALLOWED_PROGRAM_METHODS:
            raise ValueError(f"Unsupported methylation program method: {method}")
        if method in seen:
            continue
        out.append(method)
        seen.add(method)
    return out


def default_link_methods_for_preset(preset: str | None) -> tuple[str, ...]:
    normalized = normalize_program_preset(preset)
    if normalized == PROGRAM_PRESET_CONNECTABLE:
        return ("promoter_overlap", "distance_decay")
    if normalized == PROGRAM_PRESET_QC:
        return ("promoter_overlap", "nearest_tss")
    if normalized == PROGRAM_PRESET_EXPERIMENTAL:
        return ("promoter_overlap", "nearest_tss", "distance_decay")
    return ("promoter_overlap", "nearest_tss", "distance_decay")


def resolve_link_methods(link_method_expr: str | None, program_preset: str | None) -> list[str]:
    text = "all" if link_method_expr is None else str(link_method_expr).strip()
    if not text:
        text = "all"

    if text == "all":
        return list(default_link_methods_for_preset(program_preset))

    tokens = [tok.strip() for tok in text.split(",") if tok.strip()]
    if not tokens:
        tokens = ["all"]

    out: list[str] = []
    seen: set[str] = set()
    for token in tokens:
        methods = list(_LINK_METHODS) if token == "all" else [token]
        for method in methods:
            if method not in _LINK_METHODS and method != _EXTERNAL_LINK_METHOD:
                raise ValueError(f"Unsupported link_method token: {method}")
            if method in seen:
                continue
            out.append(method)
            seen.add(method)
    return out


def preferred_link_method_for_program(preset: str | None, program_method: str) -> str | None:
    normalized = normalize_program_preset(preset)
    if normalized == PROGRAM_PRESET_CONNECTABLE:
        if program_method == PROGRAM_PROMOTER_ACTIVITY:
            return "promoter_overlap"
        if program_method == PROGRAM_DISTAL_ACTIVITY:
            return "distance_decay"
        return None
    if normalized == PROGRAM_PRESET_QC:
        if program_method == PROGRAM_PROMOTER_ACTIVITY:
            return "promoter_overlap"
        if program_method == PROGRAM_LINKED_ACTIVITY:
            return "nearest_tss"
        return None
    return None


def link_methods_for_program(
    preset: str | None,
    program_method: str,
    link_methods: list[str],
) -> list[str]:
    preferred = preferred_link_method_for_program(preset, program_method)
    if preferred is None:
        return list(link_methods)
    return [method for method in link_methods if method == preferred]

