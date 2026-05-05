"""Regression tests for update14 snapshot channel/WCS exactness.

These tests intentionally avoid astropy/necstdb and check only deterministic
frequency-axis bookkeeping:
- blank optional beam fields in an embedded snapshot must not prevent snapshot
  bundle loading;
- saved_ch_start/stop must propagate to local explicit-WCS without an off-by-one
  channel shift.
"""

from types import SimpleNamespace

from tools.necst import beam_model
from tools.necst.multibeam_beam_measurement import config_separation_model as cfg


def test_snapshot_beam_blank_optional_fields_are_missing_not_numeric():
    snapshot = {
        "beams": {
            "B00": {
                "beam_id": "B00",
                "model": "legacy",
                "rotation_mode": "none",
                "az_offset_arcsec": 0,
                "el_offset_arcsec": 0,
                "reference_angle_deg": 0,
                "rotation_sign": 1,
                "rotation_slope_deg_per_deg": "",
                "dewar_angle_deg": 0,
                "pure_rotation_offset_x_el0_arcsec": "",
                "pure_rotation_offset_y_el0_arcsec": "",
                "pure_rotation_sign": "",
            }
        }
    }
    doc = beam_model.beam_model_from_snapshot(snapshot)
    b00 = doc.beams["B00"]
    assert b00["rotation_slope_deg_per_deg"] is None
    assert b00["pure_rotation_offset_x_el0_arcsec"] is None
    assert b00["pure_rotation_offset_y_el0_arcsec"] is None
    assert b00["pure_rotation_sign"] is None


def test_snapshot_saved_window_materializes_local_explicit_wcs_exactly():
    stream_id = "xffts_board4__13CO_J1_0"
    stream_entry = {
        "stream_id": stream_id,
        "frequency_axis_id": "xffts_2p5GHz_32768ch",
        "lo_chain": "rx100_13co_c18o10_usb_lsb",
        "rest_frequency_hz": 110_201_000_000.0,
    }
    layout = SimpleNamespace(saved_ch_start=7909, saved_ch_stop=9450, saved_nchan=1541, full_nchan=32768)
    snapshot = {
        "frequency_axes": {
            "xffts_2p5GHz_32768ch": {
                "if_freq_at_full_ch0_hz": 0.0,
                "if_freq_step_hz": 76296.273689992988,
                "full_nchan": 32768,
                "ctype1": "FREQ",
                "cunit1": "Hz",
                "specsys": "TOPOCENT",
                "veldef": "RADIO",
                "store_freq_column": False,
            }
        },
        "lo_chains": {
            "rx100_13co_c18o10_usb_lsb": {
                "signed_lo_sum_hz": 110_850_000_000.0,
                "if_frequency_sign": -1,
            }
        },
    }

    out = cfg._snapshot_axis_to_legacy_frequency_axis(
        stream_id=stream_id,
        stream_entry=stream_entry,
        layout=layout,
        snapshot=snapshot,
    )

    step = 76296.273689992988
    expected_first = 110_850_000_000.0 - 7909 * step
    expected_last = 110_850_000_000.0 - (9450 - 1) * step

    assert out["definition_mode"] == "explicit_wcs"
    assert out["nchan"] == 1541
    assert out["crpix1"] == 1.0
    assert out["cdelt1_hz"] == -step
    assert abs(out["crval1_hz"] - expected_first) < 1.0e-6
    local_last = out["crval1_hz"] + (out["nchan"] - 1) * out["cdelt1_hz"]
    assert abs(local_last - expected_last) < 1.0e-3
    assert out["snapshot_saved_ch_start"] == 7909
    assert out["snapshot_saved_ch_stop"] == 9450
    assert out["restfreq_hz"] == 110_201_000_000.0
