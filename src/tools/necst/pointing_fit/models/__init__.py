from __future__ import annotations

import importlib
from typing import cast

from .base import PointingModelBase
from .nanten2 import Nanten2ActualV1

_ALIASES = {
    "nanten2": Nanten2ActualV1.name,
}


def canonical_model_name(name: str) -> str:
    base = str(name or "").strip()
    if not base:
        return base
    parts = base.split("__", 1)
    head = _ALIASES.get(parts[0], parts[0])
    return head if len(parts) == 1 else f"{head}__{parts[1]}"
from .omu1p85m import (
    Omu1p85mActualV1,
    Omu1p85mCombinedReducedV1,
    Omu1p85mCombinedV1,
    Omu1p85mOpticalBasic,
    Omu1p85mRadioBasic,
)

_BUILTINS = {
    Omu1p85mActualV1.name: Omu1p85mActualV1,
    Omu1p85mCombinedV1.name: Omu1p85mCombinedV1,
    Omu1p85mCombinedReducedV1.name: Omu1p85mCombinedReducedV1,
    Omu1p85mOpticalBasic.name: Omu1p85mOpticalBasic,
    Omu1p85mRadioBasic.name: Omu1p85mRadioBasic,
    Nanten2ActualV1.name: Nanten2ActualV1,
}


def resolve_model(name: str) -> PointingModelBase:
    name = canonical_model_name(name)
    if name in _BUILTINS:
        return _BUILTINS[name]()
    if ":" in name:
        mod_name, obj_name = name.split(":", 1)
        mod = importlib.import_module(mod_name)
        cls = getattr(mod, obj_name)
        return cast(PointingModelBase, cls())
    raise KeyError(f"Unknown model: {name}")


__all__ = [
    "PointingModelBase",
    "resolve_model",
    "Omu1p85mActualV1",
    "Omu1p85mCombinedV1",
    "Omu1p85mCombinedReducedV1",
    "Omu1p85mOpticalBasic",
    "Omu1p85mRadioBasic",
    "Nanten2ActualV1",
    "canonical_model_name",
]
