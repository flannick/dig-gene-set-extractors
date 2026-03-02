"""Compatibility namespace for dig-gene-set-extractors.

This package is a neutral alias over the legacy ``omics2geneset`` package.
Use ``geneset_extractors`` for new integrations; existing ``omics2geneset``
imports remain supported.
"""

from __future__ import annotations

import importlib
import pkgutil
import sys

_LEGACY_ROOT = "omics2geneset"
_ALIASED_SUBMODULES = (
    "converters",
    "core",
    "hashing",
    "io",
    "methylation",
    "registry",
    "resource_manager",
    "rnaseq",
    "workflows",
)


def _alias_module(suffix: str) -> None:
    legacy_name = f"{_LEGACY_ROOT}.{suffix}"
    alias_name = f"{__name__}.{suffix}"
    module = importlib.import_module(legacy_name)
    sys.modules[alias_name] = module


_legacy_root_module = importlib.import_module(_LEGACY_ROOT)
__version__ = getattr(_legacy_root_module, "__version__", "0.0.0")

for _suffix in _ALIASED_SUBMODULES:
    _alias_module(_suffix)

# Alias nested modules so imports like geneset_extractors.core.validate resolve
# to the same module objects as omics2geneset.core.validate.
if hasattr(_legacy_root_module, "__path__"):
    for _module_info in pkgutil.walk_packages(_legacy_root_module.__path__, prefix=f"{_LEGACY_ROOT}."):
        _legacy_name = _module_info.name
        if _legacy_name == f"{_LEGACY_ROOT}.cli":
            continue
        _alias_name = _legacy_name.replace(f"{_LEGACY_ROOT}.", f"{__name__}.", 1)
        _module = importlib.import_module(_legacy_name)
        sys.modules[_alias_name] = _module


def __getattr__(name: str):
    return getattr(_legacy_root_module, name)


def __dir__() -> list[str]:
    return sorted(set(dir(_legacy_root_module)))
