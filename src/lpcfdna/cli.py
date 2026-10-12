"""Command line interface.

Four subcommands, in the order you would run them::

    lpcfdna sweep      # the detectability study, writes results/*.csv
    lpcfdna figures    # charts from those CSVs
    lpcfdna junction   # the copy neutral evidence calculation
    lpcfdna summary    # the headline numbers, printed

Everything is deterministic given the seeds in :class:`~lpcfdna.sweep.SweepConfig`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from .analysis import copy_neutral_summary, depth_improvement, lod_matrix, lod_table
from .figures import make_all
from .junction import CapturePanel, JunctionModel, report
from .sweep import SweepConfig, run_sweep, write_csv

REPO_ROOT = Path(__file__).resolve().parents[2]


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--results", type=Path,
                        default=REPO_ROOT / "results",
                        help="directory holding sweep.csv and friends")


def cmd_sweep(args: argparse.Namespace) -> int:
    config = SweepConfig()
    if args.quick:
        config.depths = (0.01, 0.1, 1.0)
        config.n_calibration = 200
        config.n_fpr_check = 100
        config.n_replicates = 40
        print("quick mode: fewer depths and replicates, numbers are indicative")
    print(f"{config.n_conditions} conditions, {config.n_calls} calls planned")

    rows, calibration = run_sweep(config, verbose=not args.quiet)
    write_csv(rows, args.results / "sweep.csv")
    write_csv(calibration, args.results / "calibration.csv")

    sweep = pd.DataFrame(rows)
    table = lod_table(sweep)
    table.to_csv(args.results / "limit_of_detection.csv", index=False)
    print(f"wrote {args.results / 'sweep.csv'} and two companions")
    return 0


def cmd_figures(args: argparse.Namespace) -> int:
    sweep = pd.read_csv(args.results / "sweep.csv")
    paths = make_all(sweep, args.figures)
    for path in paths:
        print(f"wrote {path}")
    return 0


def cmd_junction(args: argparse.Namespace) -> int:
    model = JunctionModel(min_reads=args.min_reads)
    panel = CapturePanel(n_fusions=args.fusions)
    text = report(model, panel)
    print(text)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n")
        print(f"\nwrote {args.out}")
    return 0


def cmd_summary(args: argparse.Namespace) -> int:
    sweep = pd.read_csv(args.results / "sweep.csv")
    calibration = pd.read_csv(args.results / "calibration.csv")

    print("CALIBRATION")
    print(calibration[["depth", "threshold", "target_fpr",
                       "achieved_fpr"]].to_string(index=False))

    print("\nLIMIT OF DETECTION (tumour fraction at 90% power)")
    print(lod_matrix(sweep).round(3).to_string())

    print("\nRELATIVE IMPROVEMENT IN DETECTION LIMIT (percent)")
    depths = sorted(sweep["depth"].unique())
    if len(depths) >= 4:
        low_pair = (depths[0], depths[2])
        high_pair = (depths[-3], depths[-1])
        print(f"  {low_pair[0]}x to {low_pair[1]}x:")
        print(depth_improvement(sweep, *low_pair).to_string())
        print(f"  {high_pair[0]}x to {high_pair[1]}x:")
        print(depth_improvement(sweep, *high_pair).to_string())

    print("\nCOPY NEUTRAL CONTROL")
    for key, value in copy_neutral_summary(sweep).items():
        print(f"  {key}: {value}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="lpcfdna",
        description="Detectability of cytogenetic lesions in low pass cell "
                    "free DNA.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_sweep = sub.add_parser("sweep", help="run the detectability sweep")
    _add_common(p_sweep)
    p_sweep.add_argument("--quick", action="store_true",
                         help="small run for smoke testing, minutes not hours")
    p_sweep.add_argument("--quiet", action="store_true")
    p_sweep.set_defaults(func=cmd_sweep)

    p_fig = sub.add_parser("figures", help="draw the figures")
    _add_common(p_fig)
    p_fig.add_argument("--figures", type=Path, default=REPO_ROOT / "figures")
    p_fig.set_defaults(func=cmd_figures)

    p_junc = sub.add_parser("junction", help="copy neutral evidence costs")
    p_junc.add_argument("--min-reads", type=int, default=3)
    p_junc.add_argument("--fusions", type=int, default=10)
    p_junc.add_argument("--out", type=Path, default=None)
    p_junc.set_defaults(func=cmd_junction)

    p_sum = sub.add_parser("summary", help="print the headline numbers")
    _add_common(p_sum)
    p_sum.set_defaults(func=cmd_summary)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
