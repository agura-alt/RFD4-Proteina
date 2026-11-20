"""Report the iPAE distribution of one or more runs, and compare arms side by side.

The campaign gates on `oracle_ppi.min_ipae_coi <= 1.2` (pass criterion `min_ipae`).
Pass/fail alone hides the thing that actually matters early on: the 2026-11-18 smoke
run failed the gate 12/12, but the useful number was that the values were 15.8-24.5 --
an order of magnitude away, not near-misses. A pass rate of 0 looks the same either way.
This prints the distribution so that difference is visible.

Usage, from the repository root:

    python CSS/report_ipae.py <out_dir> [<out_dir> ...] [--column COL] [--gate 1.2]
                              [--top N] [--csv OUT.csv]

Each <out_dir> is a run directory (the one holding inference_metrics/), and is labelled
by its directory name. Passing several compares arms in one table.

Notes on the columns:
  * `oracle_ppi.min_ipae_coi` is ONE np.min over every chain pair involving the chain of
    interest, so the tightest pair hides the rest (oracle_metric.py:136). Per-pair
    `oracle_ppi.min_ipae_<chain_i>_<chain_j>` columns are printed alongside it for that
    reason -- use --column to switch the headline to one of them.
  * `oracle_ppi.ipae_coi` is the MEAN over those pairs, not the min. Different number;
    the gate is on the min.
  * The aggregate file's `frac_under_thresholds` at 1 / 1.2 / 1.5 / 2 is recomputed here
    from the per-design rows, so it stays right after any filtering.

Design quality is not iPAE alone. A design can score well here and still fail
hotspot_adherence or the apo/holo criterion; see CSS/audit_smoke_metrics.py first to
confirm these numbers are real measurements rather than an indeterminate gate.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import pandas as pd

DEFAULT_COLUMN = "oracle_ppi.min_ipae_coi"
DEFAULT_GATE = 1.2
THRESHOLDS = [1.0, 1.2, 1.5, 2.0, 5.0, 10.0]


def load_run(out_dir: Path) -> pd.DataFrame:
    """Read every per-example metrics CSV under a run directory, newest row per design."""
    csvs = sorted(out_dir.rglob("per_example_metrics*.csv"), key=lambda f: f.stat().st_mtime)
    if not csvs:
        raise FileNotFoundError(
            f"no per_example_metrics*.csv under {out_dir}\n"
            "       Check the path is the run's out_dir (the directory containing\n"
            "       inference_metrics/), not the repository root or inference_samples/."
        )
    frames = [pd.read_csv(c) for c in csvs]
    df = pd.concat(frames, ignore_index=True)
    # Oldest->newest above, so a rerun's fresh rows win, matching aggregate_metrics.py.
    if "example_id" in df.columns:
        df = df.drop_duplicates(subset="example_id", keep="last")
    df["run"] = out_dir.name
    return df


def histogram(values: pd.Series, bins: int = 10, width: int = 40) -> list[str]:
    values = values.dropna()
    lo, hi = float(values.min()), float(values.max())
    if math.isclose(lo, hi):
        return [f"  all {len(values)} value(s) at {lo:.2f}"]
    edges = [lo + (hi - lo) * i / bins for i in range(bins + 1)]
    counts = pd.cut(values, bins=edges, include_lowest=True).value_counts().sort_index()
    peak = int(counts.max()) or 1
    lines = []
    for interval, n in counts.items():
        bar = "#" * int(round(width * n / peak))
        lines.append(f"  {interval.left:7.2f} - {interval.right:7.2f} | {bar:<{width}} {n}")
    return lines


def describe(df: pd.DataFrame, column: str, gate: float, top: int) -> bool:
    """Print one run's distribution. Returns False if the column is unusable."""
    label = df["run"].iloc[0]
    print(f"=== {label} ===")
    print(f"  {len(df)} design row(s)")
    if column not in df.columns:
        near = [c for c in df.columns if "ipae" in c.lower()]
        print(f"  MISSING column `{column}`.")
        print(f"  iPAE-ish columns present: {', '.join(near) if near else '(none)'}")
        print("  No iPAE columns at all means the oracle_ppi refold did not run or did not")
        print("  persist; CSS/audit_smoke_metrics.py reports why.\n")
        return False

    values = pd.to_numeric(df[column], errors="coerce")
    n_missing = int(values.isna().sum())
    values = values.dropna()
    if values.empty:
        print(f"  `{column}` exists but every value is NaN -- the refold produced no interface.\n")
        return False

    q = values.quantile([0.0, 0.25, 0.5, 0.75, 1.0])
    print(f"  {column}  (lower is better; gate <= {gate})")
    print(f"    n={len(values)}  missing={n_missing}  mean={values.mean():.2f}  sd={values.std():.2f}")
    print(f"    min={q[0.0]:.2f}  q25={q[0.25]:.2f}  median={q[0.5]:.2f}  q75={q[0.75]:.2f}  max={q[1.0]:.2f}")
    print()
    for line in histogram(values):
        print(line)
    print()
    print("  fraction at or under each threshold:")
    for t in THRESHOLDS:
        n = int((values <= t).sum())
        marker = "   <-- campaign gate" if math.isclose(t, gate) else ""
        print(f"    <= {t:5.2f}   {n:4d}/{len(values)}  ({100 * n / len(values):5.1f}%){marker}")

    # Per-chain-pair columns: min_ipae_coi is a single min, so a good headline can be one
    # tight pair while another interface is poor.
    pair_cols = sorted(
        c for c in df.columns
        if c.startswith("oracle_ppi.min_ipae_") and c != "oracle_ppi.min_ipae_coi"
    )
    if pair_cols:
        print("\n  per chain pair (median over designs):")
        for c in pair_cols:
            v = pd.to_numeric(df[c], errors="coerce").dropna()
            if not v.empty:
                print(f"    {c.split('.', 1)[1]:28s} median={v.median():7.2f}  min={v.min():7.2f}")

    if top:
        id_col = "example_id" if "example_id" in df.columns else df.columns[0]
        cols = [id_col, column]
        for extra in ("hotspot_obedience.hotspot_obedience", "min_ipae.pass",
                      "hotspot_adherence.pass", "ppi.pass", "antibody.pass"):
            if extra in df.columns:
                cols.append(extra)
        best = df.assign(**{column: pd.to_numeric(df[column], errors="coerce")})
        best = best.dropna(subset=[column]).nsmallest(top, column)
        print(f"\n  best {len(best)} design(s) by {column}:")
        print(best[cols].to_string(index=False, max_colwidth=44))
    print()
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dirs", nargs="+", type=Path)
    ap.add_argument("--column", default=DEFAULT_COLUMN, help=f"default {DEFAULT_COLUMN}")
    ap.add_argument("--gate", type=float, default=DEFAULT_GATE)
    ap.add_argument("--top", type=int, default=5, help="show this many best designs; 0 to skip")
    ap.add_argument("--csv", type=Path, help="write the pooled rows (id + iPAE + pass columns) here")
    args = ap.parse_args()

    frames, usable = [], 0
    for out_dir in args.out_dirs:
        try:
            df = load_run(out_dir)
        except FileNotFoundError as exc:
            print(f"=== {out_dir} ===\n  FAIL: {exc}\n")
            continue
        frames.append(df)
        if describe(df, args.column, args.gate, args.top):
            usable += 1

    if not frames:
        print("No runs could be read.")
        return 1

    if len(frames) > 1 and args.column in pd.concat(frames, ignore_index=True).columns:
        pooled = pd.concat(frames, ignore_index=True)
        pooled[args.column] = pd.to_numeric(pooled[args.column], errors="coerce")
        print("=" * 72)
        print(f"side by side on {args.column}")
        summary = pooled.groupby("run")[args.column].agg(
            n="count", min="min", median="median", max="max",
        )
        summary[f"<= {args.gate}"] = pooled.assign(
            ok=pooled[args.column] <= args.gate
        ).groupby("run")["ok"].sum()
        print(summary.to_string())
        print()
        print("Pass rates are NOT comparable across the hotspot arms: hotspot_obedience is a")
        print("fraction of conditioned hotspots, and the arms condition on 3 (fixed), 5-6")
        print("(add) and 2-3 (swap). iPAE itself is comparable; the hotspot gate is not.")

    if args.csv:
        pooled = pd.concat(frames, ignore_index=True)
        keep = ["run"] + [c for c in pooled.columns
                          if c == "example_id" or "ipae" in c.lower() or c.endswith(".pass")]
        pooled[keep].to_csv(args.csv, index=False)
        print(f"\nwrote {args.csv}")

    return 0 if usable else 1


if __name__ == "__main__":
    sys.exit(main())
