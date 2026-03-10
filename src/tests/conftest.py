from __future__ import annotations

import sys
import types
from pathlib import Path


def _add_repo_src_to_path() -> None:
    here = Path(__file__).resolve()
    repo_root = here.parents[1]
    src = repo_root / "src"
    for candidate in (src, repo_root):
        if candidate.exists():
            s = str(candidate)
            if s not in sys.path:
                sys.path.insert(0, s)


def _install_spectral_cube_stub_if_missing() -> None:
    try:
        import spectral_cube  # noqa: F401
        return
    except Exception:
        pass

    mod = types.ModuleType("spectral_cube")

    class SpectralCube:  # pragma: no cover - import stub only
        @classmethod
        def read(cls, *args, **kwargs):
            raise RuntimeError("spectral_cube stub in use; install spectral-cube for full tests")

    mod.SpectralCube = SpectralCube
    sys.modules["spectral_cube"] = mod


def _install_utils_stub_if_missing() -> None:
    try:
        import sd_radio_spectral_fits.utils  # noqa: F401
        return
    except Exception:
        pass

    mod = types.ModuleType("sd_radio_spectral_fits.utils")

    def parse_windows(windows):
        return windows

    def in_any_windows(v_axis, windows):
        import numpy as np
        mask = np.zeros(len(v_axis), dtype=bool)
        for lo, hi in windows:
            lo_f = float(lo)
            hi_f = float(hi)
            lo_f, hi_f = min(lo_f, hi_f), max(lo_f, hi_f)
            mask |= (v_axis >= lo_f) & (v_axis <= hi_f)
        return mask

    mod.parse_windows = parse_windows
    mod.in_any_windows = in_any_windows
    sys.modules["sd_radio_spectral_fits.utils"] = mod


_add_repo_src_to_path()
_install_spectral_cube_stub_if_missing()
_install_utils_stub_if_missing()
