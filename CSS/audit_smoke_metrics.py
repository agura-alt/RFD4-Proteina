"""Audit a smoke run's per-example metrics for COLUMN COVERAGE, not design quality.

This campaign's historical failure mode is not crashing jobs -- it is gates that score
every design indeterminate because the column behind them was never emitted. The
2026-11-04 run exited 0 with hotspot_adherence.all_metrics_present = 0 on 12/12 designs.
Exit status told you nothing; this script is what would have caught it.

Usage:
    python CSS/audit_smoke_metrics.py <out_dir>

Checks, for every design in per_example_metrics*.csv:
  * each required metric family emitted at least one column
  * each pass criterion reports all_metrics_present == 1 (i.e. it is really pass/fail,
    not indeterminate)
  * the three columns specifically at risk in this config are populated, not all-NaN
Exits non-zero if anything is missing.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# (label, column prefix, why it must be there)
REQUIRED_FAMILIES = [
    ("hotspot_obedience", "hotspot_obedience.", "the gate this campaign filters on"),
    ("interchain_interactions", "interchain_interactions.", "tag `interchain_interactions`"),
    ("oracle_ppi", "oracle_ppi.", "holo refold; tag `ppi_refolding`/`antibody`"),
    ("oracle_apo", "oracle_apo.apo.", "apo refold; tag `apo_refolding`"),
    ("refolded_backbone_rmsd", "refolded_backbone_rmsd.", "re-masked onto oracle_ppi"),
    ("refolded_target_aligned_binder_rmsd", "refolded_target_aligned_binder_rmsd.", "feeds antibody.pass"),
    ("refolded_target_backbone_aligned_binder_rmsd", "refolded_target_backbone_aligned_binder_rmsd.", "feeds ppi.pass"),
    ("refolded_apo_holo_rmsd", "refolded_apo_holo_rmsd.", "feeds apo_holo_consistency"),
    ("refolded_apo_holo_cdr_rmsd", "refolded_apo_holo_cdr_rmsd.", "AT RISK: align mask uses sa_mask('C_SEQ')"),
    ("antibody_contacts", "antibody_contacts.", "NEW: unlocked by the `antibody` tag"),
    ("design_motif_allatom_rmsd", "design_motif_allatom_rmsd.", "NEW: unlocked by the coordinate-conditioning tags"),
    ("sequence_adherence", "sequence_adherence.", "tag `design`"),
]

# Pass criteria that must be genuinely pass/fail rather than indeterminate.
REQUIRED_CRITERIA = [
    ("min_ipae", "this campaign's ipae gate"),
    ("hotspot_adherence", "all three hotspots contacted"),
    ("apo_holo_consistency", "apo/holo self-consistency"),
    ("ppi", "repo preset, needs the re-masked backbone-aligned RMSD"),
    ("antibody", "NEW: unlocked by the `antibody` tag; depends on an at-risk mask"),
]

# Columns that are expected to carry real numbers on at least one design.
MUST_BE_POPULATED = [
    "refolded_apo_holo_cdr_rmsd",
    "design_motif_allatom_rmsd",
    "antibody_contacts",
]


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    out_dir = Path(sys.argv[1])
    csvs = sorted(out_dir.rglob("per_example_metrics*.csv"))
    if not csvs:
        print(f"FAIL: no per_example_metrics*.csv under {out_dir}")
        return 1

    df = pd.concat([pd.read_csv(c) for c in csvs], ignore_index=True)
    print(f"Read {len(csvs)} file(s), {len(df)} design row(s), {len(df.columns)} columns")
    print(f"  from: {out_dir}\n")

    problems: list[str] = []

    print("--- metric families ---")
    for label, prefix, why in REQUIRED_FAMILIES:
        cols = [c for c in df.columns if c.startswith(prefix)]
        if cols:
            print(f"  ok      {label:46s} {len(cols)} column(s)")
        else:
            print(f"  MISSING {label:46s} ({why})")
            problems.append(f"no columns for {label} -- {why}")

    print("\n--- pass criteria (must be determinate) ---")
    for label, why in REQUIRED_CRITERIA:
        pass_col, present_col = f"{label}.pass", f"{label}.all_metrics_present"
        if pass_col not in df.columns:
            print(f"  MISSING {label:46s} no `{pass_col}` column ({why})")
            problems.append(f"no {pass_col} -- {why}")
            continue
        if present_col in df.columns:
            frac = float(pd.to_numeric(df[present_col], errors="coerce").fillna(0).mean())
            if frac < 1.0:
                print(f"  INDET.  {label:46s} all_metrics_present={frac:.2f} on {len(df)} design(s)")
                problems.append(f"{label} is indeterminate (all_metrics_present={frac:.2f})")
                continue
        n_pass = int(pd.to_numeric(df[pass_col], errors="coerce").fillna(0).sum())
        print(f"  ok      {label:46s} determinate; {n_pass}/{len(df)} passed")

    print("\n--- at-risk columns must carry real values ---")
    for prefix in MUST_BE_POPULATED:
        cols = [c for c in df.columns if c.startswith(prefix + ".")]
        if not cols:
            continue  # already reported above
        numeric = df[cols].apply(pd.to_numeric, errors="coerce")
        if numeric.notna().any().any():
            print(f"  ok      {prefix:46s} has non-NaN values")
        else:
            print(f"  ALL-NaN {prefix:46s} columns exist but every value is NaN")
            problems.append(f"{prefix} columns are entirely NaN")

    print("\n" + "=" * 70)
    if problems:
        print(f"{len(problems)} PROBLEM(S) -- do not scale this campaign yet:")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("Metric coverage is complete. Design quality is a separate question;")
    print("a smoke run of this size says nothing about pass rates.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
