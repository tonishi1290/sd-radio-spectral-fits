# Multi-beam beam measurement package

This package is an additive implementation around the existing `sun_scan_v4_v5.py`
and `necst_v4_sdfits_converter.py` scripts.

## Included files

- `sun_scan_v4_v5.py` — original script, kept unchanged
- `necst_v4_sdfits_converter.py` — original converter, kept unchanged
- `multibeam_beam_measurement/` — reusable package
- `sunscan_extract_multibeam.py` — convenience CLI wrapper
- `sunscan_fit_multibeam.py` — convenience CLI wrapper
- `sunscan_multibeam.py` — subcommand wrapper

## Basic workflow

1. Verify that the original `sun_scan_v4_v5.py` still works in your environment.
2. Verify that the original `necst_v4_sdfits_converter.py` still works in your environment.
3. Run multi-beam extraction:

```bash
python sunscan_extract_multibeam.py RAWDATA \
  --spectrometer-config beams.toml \
  --outdir out_extract
```

4. Run beam fitting:

```bash
python sunscan_fit_multibeam.py \
  out_extract/sunscan_multibeam_scan_summary_<TAG>.csv \
  --spectrometer-config beams.toml \
  --outdir out_fit \
  --center-beam-id B00
```

## Notes

- Per-stream `sun_scan`-compatible outputs are written under `outdir/per_stream/<stream_name>/`.
- Aggregated extraction output is written as `sunscan_multibeam_scan_summary_<tag>.csv`.
- Fit output includes `beam_fit_results_<model>.csv`, `run_shift_results_<model>.csv`,
  `beam_model_<model>.toml`, and `fit_summary.txt`.
- The implementation intentionally keeps the original scripts untouched to reduce risk.
