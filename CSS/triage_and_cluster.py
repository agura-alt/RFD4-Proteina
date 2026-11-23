#!/usr/bin/env python
"""Triage the confident-interface designs and cluster them by structure.

Run on the cluster that holds the output trees: `designed_chain_path` in the
metrics CSVs is an absolute path on that machine, so clustering elsewhere would
need the column rewritten first.

    python triage_and_cluster.py <run_root>/*/inference_metrics/per_example_metrics.csv

Writes hits_clustered.csv: one row per design at or under the 5 A
confident-interface cut, with a cluster id at three TM-score thresholds.

To pool an arm whose outputs came from another cluster, copy its chain CIFs onto
this machine and rebase the paths recorded in its CSV:

    --remap /anvil/scratch/.../runs=/delta/imported/anvil_runs

Every designed_chain_path starting with the left side is rewritten to the right.
Repeat the flag for more than one source. Paths that do not exist after remapping
are reported and dropped rather than silently clustered as absent.

Needs FOLDSEEK_EXEC set and the chain-of-interest CIFs on disk (generation with
save_chain_of_interest=true, and aggregation not having cleaned them up).
"""

import sys
from pathlib import Path

import pandas as pd

from rfproteina.inference_utils.foldseek import diversity_foldseek

IPAE = "oracle_ppi.min_ipae_coi"
CUT = 5.0
# 0.5 is the pipeline default and asks "same fold?" -- too coarse for a panel of
# VHHs on one scaffold. 0.7/0.8 are what separate CDR-loop solutions.
THRESHOLDS = (0.5, 0.7, 0.8)


def main(argv: list[str]) -> int:
    remaps: list[tuple[str, str]] = []
    paths: list[str] = []
    i = 0
    while i < len(argv):
        if argv[i] == "--remap":
            old, _, new = argv[i + 1].partition("=")
            if not new:
                print("--remap takes OLD_PREFIX=NEW_PREFIX")
                return 2
            remaps.append((old, new))
            i += 2
        else:
            paths.append(argv[i])
            i += 1

    if not paths:
        print(__doc__)
        return 2

    frames = []
    for p in paths:
        d = pd.read_csv(p)
        # <arm>/inference_metrics/per_example_metrics.csv -> <arm>
        d["arm"] = Path(p).parent.parent.name
        frames.append(d)
    df = pd.concat(frames, ignore_index=True)

    if IPAE not in df.columns:
        print(f"No {IPAE} column. Columns containing 'ipae': "
              f"{[c for c in df.columns if 'ipae' in c]}")
        return 1

    print(f"{len(df)} designs, {int(df[IPAE].notna().sum())} with an iPAE measurement")
    hits = df[df[IPAE] <= CUT].copy()
    print(f"{len(hits)} at or under {CUT} A")
    print(hits.groupby("arm").size().to_string())
    if hits.empty:
        return 1

    if "designed_chain_path" not in hits.columns:
        print("No designed_chain_path column -- the runs were generated without "
              "save_chain_of_interest=true. Re-run rfd4-evaluate over the "
              "complexes to write the chain CIFs, then re-run this.")
        return 1

    missing = hits["designed_chain_path"].isna().sum()
    if missing:
        print(f"WARNING: dropping {missing} hit(s) with no chain CIF from clustering")
    hits = hits[hits["designed_chain_path"].notna()].copy()

    if remaps:
        def rebase(p: str) -> str:
            for old, new in remaps:
                if p.startswith(old):
                    return new + p[len(old):]
            return p
        hits["designed_chain_path"] = hits["designed_chain_path"].map(rebase)

    # A path that is merely absent would otherwise reach foldseek as a missing
    # file and take its design out of the clustering without comment.
    exists = hits["designed_chain_path"].map(lambda p: Path(p).is_file())
    if not exists.all():
        print(f"WARNING: {int((~exists).sum())} chain CIF(s) not readable on this "
              f"machine; dropping them. First few:")
        for p in hits.loc[~exists, "designed_chain_path"].head(3):
            print(f"  {p}")
        print("  (cross-cluster rows need --remap and the files copied over)")
        hits = hits[exists]
    if hits.empty:
        return 1

    for tm in THRESHOLDS:
        ids = diversity_foldseek(hits["designed_chain_path"].tolist(), tm_threshold=tm)
        if ids is None:
            print(f"Foldseek failed at TM {tm} -- check FOLDSEEK_EXEC and that the "
                  f"chain CIFs still exist.")
            return 1
        hits[f"cluster_tm{tm}"] = ids
        print(f"TM {tm}: {len(set(ids))} cluster(s) over {len(ids)} designs")

    # Quality columns for the triage, whichever of them this run produced.
    quality = [c for c in hits.columns
               if any(t in c for t in ("ptm_coi", "plddt_coi", "clash", "hotspot"))]
    cols = (["arm", "example_id", IPAE] + quality
            + [f"cluster_tm{t}" for t in THRESHOLDS] + ["designed_chain_path"])
    cols = [c for c in cols if c in hits.columns]

    out = Path("hits_clustered.csv")
    hits.sort_values(IPAE)[cols].to_csv(out, index=False)
    print(f"Wrote {out.resolve()}")
    print("\nPick one representative per cluster at TM 0.8, ranked on apo pTM "
          "(does the binder fold alone) before iPAE.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
