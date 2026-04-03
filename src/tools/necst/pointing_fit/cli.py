from __future__ import annotations

import argparse

from .apply import ApplyOptions, run_apply
from .fit import FitOptions, run_fit
from .inspect import InspectOptions, run_inspect
from .normalize import NormalizeOptions, run_normalize
from .report import ReportOptions, run_report


def add_normalize(sub):
    p = sub.add_parser("normalize", help="normalize raw dx/dy CSV into canonical CSV")
    p.add_argument("--input", dest="input_csv")
    p.add_argument("--manifest")
    p.add_argument("--output", required=True)
    p.add_argument("--pair-mode", choices=["union", "raw"], default="union")
    p.add_argument("--dataset-id")
    p.add_argument("--used-param")
    p.add_argument("--telescope")
    p.add_argument("--model", default="omu1p85m_combined_v1")
    p.add_argument("--angle-unit", choices=["deg", "rad"], default="deg", help="fallback unit for raw Az/El columns when the column names do not already encode units (for example az_rad or el_deg)")
    p.add_argument("--offset-unit", choices=["arcsec", "deg", "rad"], default="arcsec", help="fallback unit for raw dx/dy columns when the column names do not already encode units; default is arcsec for backward compatibility with legacy d_param.csv files")
    p.add_argument("--measurement-space", choices=["model_delta", "command_offset", "star_offset"], default="star_offset", help="meaning of raw dx/dy in the input CSV; user-facing optical inputs should normally use star_offset")
    p.add_argument("--positive-az-moves-star", choices=["left", "right"], default=None, help="for star_offset inputs: where the star moves in the image when +Az command is applied; not used for command_offset/model_delta")
    p.add_argument("--positive-el-moves-star", choices=["up", "down"], default=None, help="for star_offset inputs: where the star moves in the image when +El command is applied; not used for command_offset/model_delta")
    p.add_argument("--command-err-mode", choices=["subtract", "add"], default=None, help="relation between command coordinates and model error terms: subtract => cmd=true-err, add => cmd=true+err")
    p.set_defaults(func=cmd_normalize)


def add_fit(sub):
    p = sub.add_parser("fit", help="fit telescope pointing model")
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--model", default="omu1p85m_combined_v1")
    p.add_argument("--fit-target", choices=["absolute", "delta"], default="absolute")
    p.add_argument("--solve-mode", choices=["joint", "separate"], default="joint")
    p.add_argument("--robust", dest="robust_loss", choices=["linear", "soft_l1", "huber", "cauchy"], default="soft_l1")
    p.add_argument("--f-scale", default="auto")
    p.add_argument("--clip-sigma", type=float, default=0.0)
    p.add_argument("--clip-max-iter", type=int, default=0)
    p.add_argument("--bootstrap", type=int, default=0)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument(
        "--fix",
        action="append",
        default=[],
        help="fix parameter, e.g. --fix d1[arcsec]=8.0 or --fix g1=0.0. Repeatable.",
    )
    p.add_argument(
        "--fix-file",
        default=None,
        help="TOML file containing [fixed_params] or [pointing_params] entries to hold fixed during fit.",
    )
    p.set_defaults(func=cmd_fit)


def add_apply(sub):
    p = sub.add_parser("apply", help="apply delta params to base params (delta TOML or delta-fit output directory)")
    p.add_argument("--base-param", "--base", dest="base_param", required=True)
    p.add_argument("--delta-param", "--update-param", dest="delta_param", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_apply)


def add_inspect(sub):
    p = sub.add_parser("inspect", help="inspect normalized CSV")
    p.add_argument("--input", required=True)
    p.set_defaults(func=cmd_inspect)


def add_report(sub):
    p = sub.add_parser("report", help="create plots/markdown from residuals.csv")
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--title", default="pointing-fit report")
    p.set_defaults(func=cmd_report)


def cmd_normalize(args):
    df = run_normalize(
        NormalizeOptions(
            input_csv=args.input_csv,
            manifest=args.manifest,
            output_csv=args.output,
            pair_mode=args.pair_mode,
            dataset_id=args.dataset_id,
            used_param=args.used_param,
            telescope=args.telescope,
            model=args.model,
            angle_unit=args.angle_unit,
            offset_unit=args.offset_unit,
            measurement_space=args.measurement_space,
            positive_az_moves_star=args.positive_az_moves_star,
            positive_el_moves_star=args.positive_el_moves_star,
            command_err_mode=args.command_err_mode,
        )
    )
    print(f"Wrote normalized CSV: {args.output} ({len(df)} rows)")


def cmd_fit(args):
    res = run_fit(
        FitOptions(
            input_csv=args.input,
            outdir=args.outdir,
            model=args.model,
            fit_target=args.fit_target,
            solve_mode=args.solve_mode,
            robust_loss=args.robust_loss,
            f_scale=args.f_scale,
            clip_sigma=args.clip_sigma,
            clip_max_iter=args.clip_max_iter,
            bootstrap=args.bootstrap,
            seed=args.seed,
            fixed_param_specs=list(args.fix or []),
            fixed_param_file=args.fix_file,
        )
    )
    print(f"Wrote fit outputs to: {args.outdir}")
    print(f"Residual radial RMS [arcsec]: {res['quality']['residual_rms_arcsec']['radial']:.3f}")
    if res.get("fixed_parameters"):
        print(f"Fixed parameters: {', '.join(sorted(res['fixed_parameters']))}")


def cmd_apply(args):
    out = run_apply(ApplyOptions(base_param=args.base_param, delta_param=args.delta_param, output=args.output))
    print(f"Wrote applied params: {args.output}")
    print(f"Parameters written: {len(out)}")


def cmd_inspect(args):
    print(run_inspect(InspectOptions(input_csv=args.input)))


def cmd_report(args):
    files = run_report(ReportOptions(input_csv=args.input, outdir=args.outdir, title=args.title))
    print(f"Wrote: {files['png']}")
    if files.get("fit_map_png"):
        print(f"Wrote: {files['fit_map_png']}")
    print(f"Wrote: {files['markdown']}")


def main() -> None:
    parser = argparse.ArgumentParser(prog="pointing-fit")
    sub = parser.add_subparsers(dest="command", required=True)
    add_normalize(sub)
    add_fit(sub)
    add_apply(sub)
    add_inspect(sub)
    add_report(sub)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
