from __future__ import annotations

import numpy as np
import pytest

from sd_radio_spectral_fits.map.cube_baseline.session import BaselineSession
from sd_radio_spectral_fits.map_3d.cube_baseline import orchestrator as orch



def test_session_axis_order_hint_middle_axis_standardizes_correctly() -> None:
    cube = np.arange(2 * 5 * 3, dtype=np.float32).reshape(2, 5, 3)  # (ny, nchan, nx)
    v_axis = np.arange(5, dtype=float)
    session = BaselineSession(cube, v_axis, axis_order_hint="y_v_x")
    assert session.axis_order_in == "y_v_x"
    assert session.shape == (5, 2, 3)
    # spectral axis must now be axis 0
    assert np.allclose(session.cube_original[:, 0, 0], cube[0, :, 0])



def test_session_set_prior_rejects_length_mismatch() -> None:
    session = BaselineSession(np.zeros((4, 2, 2), dtype=np.float32), np.arange(4, dtype=float))
    with pytest.raises(ValueError, match="shape mismatch"):
        session.set_prior(linefree_mask_1d=np.ones(3, dtype=bool))



def test_run_one_iteration_prior_linefree_without_auto_raises() -> None:
    session = BaselineSession(np.zeros((4, 2, 2), dtype=np.float32), np.arange(4, dtype=float))
    with pytest.raises(ValueError, match="requires automatic line-free estimation"):
        orch.run_one_iteration(
            session,
            auto_linefree=False,
            linefree_mode="prior",
            enable_ripple=False,
        )



def test_run_one_iteration_ripple_prior_falls_back_to_auto(monkeypatch) -> None:
    session = BaselineSession(np.zeros((4, 2, 2), dtype=np.float32), np.arange(4, dtype=float))

    monkeypatch.setattr(orch, "estimate_linefree_mask_from_cube", lambda cube, cfg, agg="median": np.ones(session.nchan, dtype=bool))
    monkeypatch.setattr(orch, "estimate_ripple_frequencies_fft", lambda spec, lf, rcfg, poly_order_pre: [0.125])

    called = {}

    def fake_fit_cube_baseline(cube_original, *, linefree_mask_1d, baseline_cfg, ripple_freqs, target_mask_2d):
        called["linefree"] = np.asarray(linefree_mask_1d).copy()
        called["ripple_freqs"] = list(ripple_freqs) if ripple_freqs is not None else None
        stats = {"rms_map": np.zeros((cube_original.shape[1], cube_original.shape[2]), dtype=np.float32)}
        return np.zeros_like(cube_original), stats

    monkeypatch.setattr(orch, "fit_cube_baseline", fake_fit_cube_baseline)

    orch.run_one_iteration(
        session,
        auto_linefree=True,
        linefree_mode="auto",
        enable_ripple=True,
        ripple_mode="prior",
    )

    assert np.all(called["linefree"])
    assert called["ripple_freqs"] == [0.125]



def test_manual_velocity_windows_are_excluded_from_linefree(monkeypatch) -> None:
    session = BaselineSession(np.zeros((4, 2, 2), dtype=np.float32), np.arange(4, dtype=float))

    monkeypatch.setattr(orch, "estimate_linefree_mask_from_cube", lambda cube, cfg, agg="median": np.ones(session.nchan, dtype=bool))
    monkeypatch.setattr(orch, "create_manual_signal_mask_1d", lambda v_axis, wins: np.array([False, True, False, False], dtype=bool))

    called = {}

    def fake_fit_cube_baseline(cube_original, *, linefree_mask_1d, baseline_cfg, ripple_freqs, target_mask_2d):
        called["linefree"] = np.asarray(linefree_mask_1d).copy()
        stats = {"rms_map": np.zeros((cube_original.shape[1], cube_original.shape[2]), dtype=np.float32)}
        return np.zeros_like(cube_original), stats

    monkeypatch.setattr(orch, "fit_cube_baseline", fake_fit_cube_baseline)

    orch.run_one_iteration(
        session,
        auto_linefree=True,
        linefree_mode="auto",
        manual_v_windows=[(1.0, 1.0)],
        enable_ripple=False,
    )

    assert called["linefree"].tolist() == [True, False, True, True]
