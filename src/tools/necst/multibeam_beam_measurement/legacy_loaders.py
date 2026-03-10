from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Optional

_SUN_SCAN_CACHE: Optional[ModuleType] = None
_CONVERTER_CACHE: Optional[ModuleType] = None


def _load_module_from_path(path: Path, module_name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"failed to build import spec for {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def default_sun_scan_path() -> Path:
    return Path(__file__).resolve().parent.parent / "sun_scan_v4_v5.py"


def default_converter_path() -> Path:
    return Path(__file__).resolve().parent.parent / "necst_v4_sdfits_converter.py"


def load_sun_scan_module(path: Optional[Path] = None) -> ModuleType:
    global _SUN_SCAN_CACHE
    path = Path(path) if path is not None else default_sun_scan_path()
    if _SUN_SCAN_CACHE is not None and Path(getattr(_SUN_SCAN_CACHE, "__file__", "")) == path:
        return _SUN_SCAN_CACHE
    try:
        _SUN_SCAN_CACHE = _load_module_from_path(path, "legacy_sun_scan_v4_v5")
        return _SUN_SCAN_CACHE
    except Exception as exc:
        raise RuntimeError(
            f"failed to load legacy sun_scan module from '{path}'. "
            "Please ensure the original sun_scan_v4_v5.py and its runtime dependencies are available."
        ) from exc


def load_converter_module(path: Optional[Path] = None) -> ModuleType:
    global _CONVERTER_CACHE
    path = Path(path) if path is not None else default_converter_path()
    if _CONVERTER_CACHE is not None and Path(getattr(_CONVERTER_CACHE, "__file__", "")) == path:
        return _CONVERTER_CACHE
    try:
        _CONVERTER_CACHE = _load_module_from_path(path, "legacy_necst_v4_sdfits_converter")
        return _CONVERTER_CACHE
    except Exception as exc:
        raise RuntimeError(
            f"failed to load legacy converter module from '{path}'. "
            "Please ensure necst_v4_sdfits_converter.py and its runtime dependencies are available."
        ) from exc
