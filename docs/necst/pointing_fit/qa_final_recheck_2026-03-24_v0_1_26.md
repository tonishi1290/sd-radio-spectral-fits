# pointing_fit v0.1.26 final recheck

## Additional fixes found during another deep audit

1. `normalize` used undocumented hard requirements on raw CSV columns `dt_cross_id` and `cross_id`. These are now auto-generated from row index when absent.
2. Compatibility alias module execution `python -m optical_pointing_cli.cli --help` previously produced no output because the wrapper lacked a `__main__` handoff. This now forwards to `pointing_fit.cli:main`.

## Rechecks performed

- `compileall` on `pointing_fit`, `optical_pointing_cli`, `validation`, `examples`
- local `pip install --no-deps --target ...`
- `python -m pointing_fit.cli --help` and `python -m optical_pointing_cli.cli --help`
- synthetic example `normalize -> inspect -> fit(abs/delta) -> report -> apply`
- bundled NANTEN2 manifest normalize
- custom simulations including missing `dt_cross_id`/`cross_id`, `dx_rad`/`dy_rad`, `az_rad`/`el_rad`, and sign conventions
- validation scripts

## Outcome

The package passed the above checks after the two fixes.
