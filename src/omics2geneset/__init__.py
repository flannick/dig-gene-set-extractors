"""Legacy compatibility namespace.

`omics2geneset` remains import-compatible, but canonical modules now live under
`geneset_extractors`.
"""

from __future__ import annotations

import importlib
import pkgutil
import sys

_CANON_ROOT = "geneset_extractors"
_ALIAS_MAP = {
    "registry": "registry",
    "resource_manager": "resource_manager",
    "hashing": "hashing",
    "core": "core",
    "io": "io",
    "workflows": "workflows",
    "converters": "extractors.converters",
    "atac": "extractors.atac",
    "rnaseq": "extractors.rnaseq",
    "methylation": "extractors.methylation",
    "cnv": "extractors.cnv",
    "chipseq": "extractors.chipseq",
    "proteomics": "extractors.proteomics",
}

_canon_root_module = importlib.import_module(_CANON_ROOT)
__version__ = getattr(_canon_root_module, "__version__", "0.0.0")

for legacy_suffix, canon_suffix in _ALIAS_MAP.items():
    module = importlib.import_module(f"{_CANON_ROOT}.{canon_suffix}")
    sys.modules[f"{__name__}.{legacy_suffix}"] = module


def _legacy_name_for_canonical(name: str) -> str:
    legacy_name = name.replace(f"{_CANON_ROOT}.", f"{__name__}.", 1)
    legacy_name = legacy_name.replace(f"{__name__}.extractors.converters", f"{__name__}.converters")
    legacy_name = legacy_name.replace(f"{__name__}.extractors.atac", f"{__name__}.atac")
    legacy_name = legacy_name.replace(f"{__name__}.extractors.rnaseq", f"{__name__}.rnaseq")
    legacy_name = legacy_name.replace(f"{__name__}.extractors.methylation", f"{__name__}.methylation")
    legacy_name = legacy_name.replace(f"{__name__}.extractors.cnv", f"{__name__}.cnv")
    legacy_name = legacy_name.replace(f"{__name__}.extractors.chipseq", f"{__name__}.chipseq")
    legacy_name = legacy_name.replace(f"{__name__}.extractors.proteomics", f"{__name__}.proteomics")
    return legacy_name


if hasattr(_canon_root_module, "__path__"):
    for module_info in pkgutil.walk_packages(_canon_root_module.__path__, prefix=f"{_CANON_ROOT}."):
        canon_name = module_info.name
        if canon_name == f"{_CANON_ROOT}.cli":
            continue
        module = importlib.import_module(canon_name)
        sys.modules[_legacy_name_for_canonical(canon_name)] = module


def __getattr__(name: str):
    return getattr(_canon_root_module, name)


def __dir__() -> list[str]:
    return sorted(set(dir(_canon_root_module)))
