from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import pandas as pd

from pointing_fit.io_utils import write_pointing_param_toml
from pointing_fit.models import resolve_model


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="example_out", help="output directory")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    model = resolve_model("omu1p85m_combined_v1")
    true_params = {
        "a1": 3.0 / 3600,
        "a2": -1.5 / 3600,
        "a3": 0.8 / 3600,
        "b1": 2.0 / 3600,
        "b2": -1.0 / 3600,
        "b3": 5.0 / 3600,
        "g1": 0.0,
        "c1": 1.5 / 3600,
        "c2": -2.2 / 3600,
        "d1": 8.0 / 3600,
        "d2": -6.0 / 3600,
    }
    used_params = {k: 0.0 for k in true_params}
    az = np.linspace(0, 330, 60)
    el = 25 + 45 * (0.5 + 0.5 * np.sin(np.deg2rad(az * 1.3)))
    daz = model.predict_dx_deg(true_params, az, el)
    ddel = model.predict_dy_deg(true_params, az, el)
    rng = np.random.default_rng(0)
    daz = daz + rng.normal(0, 1.0 / 3600.0, size=daz.size)
    ddel = ddel + rng.normal(0, 1.0 / 3600.0, size=ddel.size)
    df = pd.DataFrame(
        {
            "dt_cross_id": [f"c{i:03d}" for i in range(len(az))],
            "cross_id": np.arange(len(az)),
            "Az": az,
            "El": el,
            "dx_arcsec": daz * 3600.0,
            "dy_arcsec": ddel * 3600.0,
        }
    )
    df.to_csv(outdir / "d_param.csv", index=False)
    write_pointing_param_toml(outdir / "pointing_param.toml", model_name="omu1p85m_combined_v1", params_deg=used_params, unit="deg")
    manifest = '''[defaults]
angle_unit = "deg"
offset_unit = "arcsec"
measurement_space = "model_delta"

[[datasets]]
dataset_id = "synthetic"
input = "d_param.csv"
used_param = "pointing_param.toml"
telescope = "omu1p85m"
model = "omu1p85m_combined_v1"
'''
    (outdir / "datasets.toml").write_text(manifest, encoding="utf-8")
    print(f"Wrote {outdir}/d_param.csv, pointing_param.toml, datasets.toml")


if __name__ == "__main__":
    main()
