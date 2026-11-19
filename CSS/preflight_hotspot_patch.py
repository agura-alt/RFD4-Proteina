"""Offline preflight for the hotspot-sampling arms: does the resampling actually fit?

Two arms sample hotspots from a patch around the three selected atoms
(A/204 LEU CD2, A/206 LEU CD2, A/246 TYR CE1), once per replicate:

    ..._hotspotadd_generate.yaml   SampleNearbyHotspots wrapped in KeepSeedHotspots:
                                   the three seeds PLUS 2-3 sampled, so 5-6 in all.
    ..._hotspotswap_generate.yaml  SampleNearbyHotspots bare: 2-3 sampled INSTEAD of
                                   the seeds, which the sampler excludes.

Both shapes are understood here. The report says which one the config composed and
whether the seeds survived each draw, so an arm cannot quietly turn into the other.
Pass --config to choose; the default is the add arm.

That transform RAISES when it cannot fit the requested number of separated
hotspots in the patch (hotspot_sampling.py:117). On a groove epitope like this
one that is a real risk: the seeds span 13.77 A, so a 15 A-diameter sphere
centred on their centroid is not obviously roomy once seed residues are excluded,
atoms must be >= min_distance apart, and each must clear the SASA floor. A raise
there lands after the GPU is already allocated.

This script answers it without a GPU, a checkpoint or an AF3 call: it applies the
transform to each built input many times and reports the candidate pool, the
count distribution, which residues get picked, and any failures.

It instantiates the transform FROM the arm's composed config, so the draws reported
here are the draws the run will make. That doubles as the check that the
`generation.dataset.sampling_transforms` block composed at all: a mis-addressed
key is silently inert -- generate.py:204 simply finds no `sampling_transforms` and
the run proceeds with the fixed seed hotspots -- and no metric downstream would
say so.

The transform is applied to the input CIF as written. In the real pipeline it runs
after RealizeExpandableSegments, which only expands the binder's CDR sentinels; the
target chain sampled here is untouched by that step, so the target-side result is
the same.

Usage, from the repository root:
    python CSS/preflight_hotspot_patch.py [n_draws] [--config NAME]
        [--shape sphere|oval] [--patch-diameter A] [--major-diameter A]
        [--min-distance A] [--min-sasa A2] [--num-hotspots 2,3]

Flags override the composed values, for trying settings before editing the config.
Exits non-zero if the config carries no SampleNearbyHotspots, if any draw fails, or
if any input offers fewer than 3 distinct eligible atoms.
"""

from __future__ import annotations

import argparse
import collections
import os
from pathlib import Path

import hydra
import numpy as np
from hydra import compose, initialize_config_dir
from omegaconf import OmegaConf

import rfproteina.conditions  # noqa: F401  -- registers condition classes before parse
from atomworks.io import parse
from atomworks.io.utils.atom_array_plus import as_atom_array_plus
from rfproteina.conditions import C_HOT
from rfproteina.utils.parse_presets import ANNOTATED_CIF

ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "CSS" / "configs"
DEFAULT_CONFIG_NAME = "css_cdr123_L204CD2_L206CD2_Y246CE1_hotspotadd_generate"
INPUT_DIR = ROOT / "CSS" / "inputs" / "cifs_L204CD2_L206CD2_Y246CE1"
SAMPLER = "rfproteina.transforms.hotspot_sampling.SampleNearbyHotspots"
WRAPPER = "rfproteina.transforms.sampling.css_seeded_hotspots.KeepSeedHotspots"
SEED = 42  # the smoke test pins seed=42; same stream, so the same draws


def transform_entry_from_config(config_name: str) -> dict:
    """Compose the arm's config and pull out its hotspot-sampling entry, as written."""
    os.environ.setdefault("PROJECT_ROOT", str(ROOT))
    os.environ.setdefault("DATA_PATH", "/unused/for/preflight")
    with initialize_config_dir(config_dir=str(CONFIG_DIR), version_base=None):
        cfg = compose(config_name=config_name)
    node = OmegaConf.select(cfg, "generation.dataset.sampling_transforms")
    entries = OmegaConf.to_container(node, resolve=True) if node else []
    matching = [e for e in entries if isinstance(e, dict) and e.get("_target_") in (SAMPLER, WRAPPER)]
    if not matching:
        raise SystemExit(
            f"FAIL: {config_name} composes no hotspot sampler under\n"
            "      generation.dataset.sampling_transforms. The run would quietly use the\n"
            "      fixed seed hotspots instead, and nothing downstream would say so."
        )
    if len(matching) > 1:
        raise SystemExit(f"FAIL: {config_name} declares {len(matching)} hotspot-sampling transforms.")
    entry = matching[0]
    if entry["_target_"] == WRAPPER and entry.get("sampler", {}).get("_target_") != SAMPLER:
        raise SystemExit(f"FAIL: {WRAPPER} in {config_name} does not wrap a {SAMPLER}.")
    return entry


def load(path: Path):
    return as_atom_array_plus(parse(str(path), config=ANNOTATED_CIF)["assemblies"]["1"][0])


def describe(array, mask) -> list[str]:
    idx = np.flatnonzero(mask)
    return [f"{array.chain_id[i]}/{int(array.res_id[i])} {array.res_name[i]} {array.atom_name[i]}" for i in idx]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("n_draws", nargs="?", type=int, default=50)
    ap.add_argument("--config", default=DEFAULT_CONFIG_NAME)
    ap.add_argument("--shape", choices=("sphere", "oval"))
    ap.add_argument("--patch-diameter", type=float)
    ap.add_argument("--major-diameter", type=float)
    ap.add_argument("--min-distance", type=float)
    ap.add_argument("--min-sasa", type=float)
    ap.add_argument("--num-hotspots")
    args = ap.parse_args()

    entry = transform_entry_from_config(args.config)
    keeps_seeds = entry["_target_"] == WRAPPER
    inner = entry["sampler"] if keeps_seeds else entry
    print(f"config: {args.config}")
    print(f"  composed transform: {entry['_target_'].rsplit('.', 1)[-1]}"
          + (f" wrapping {SAMPLER.rsplit('.', 1)[-1]}" if keeps_seeds else ""))
    print(f"  sampler settings: {({k: v for k, v in inner.items() if k != '_target_'})}")
    print("  seeds: KEPT, sampled atoms are additional" if keeps_seeds
          else "  seeds: REPLACED by the draw (the sampler excludes the seed residues)")

    overrides = {
        "shape": args.shape,
        "patch_diameter": args.patch_diameter,
        "major_diameter": args.major_diameter,
        "min_distance": args.min_distance,
        "min_atom_sasa": args.min_sasa,
        "num_hotspots": [int(n) for n in args.num_hotspots.split(",")] if args.num_hotspots else None,
    }
    overrides = {k: v for k, v in overrides.items() if v is not None}
    if overrides:
        print(f"  COMMAND-LINE OVERRIDES (not what the run will use): {overrides}")
        inner.update(overrides)

    # Instantiate exactly as Hydra will at run time (generate.py:204), wrapper included.
    try:
        transform = hydra.utils.instantiate(OmegaConf.create(entry))
    except Exception as exc:  # noqa: BLE001 -- the message matters more than the class
        missing = entry["_target_"].rsplit(".", 1)[0]
        raise SystemExit(
            f"FAIL: could not instantiate {entry['_target_']}\n"
            f"       {type(exc).__name__}: {exc}\n"
            "\n"
            f"      If the module is {WRAPPER.rsplit('.', 1)[0]}, check it exists in THIS\n"
            "      checkout: src/rfproteina/transforms/sampling/ is git-ignored (.gitignore:9),\n"
            "      so a checkout made by `git clone` or `git pull` does not have it. Copy the\n"
            "      module across, or force-add it so it travels with the config that names it:\n"
            f"          git add -f {missing.replace('.', '/').replace('rfproteina/', 'src/rfproteina/', 1)}.py"
        ) from exc
    sampler = transform.sampler if keeps_seeds else transform
    # Read the radius back off the instance, so anything the config left at the
    # transform's own default is still reported correctly below.
    radius = (sampler.major_diameter if sampler.shape == "oval" else sampler.patch_diameter) / 2
    print(f"\n{args.n_draws} draws per input, numpy seed {SEED}\n")

    cifs = sorted(INPUT_DIR.glob("*.cif"))
    if not cifs:
        print(f"FAIL: no .cif inputs under {INPUT_DIR}")
        return 1

    problems: list[str] = []
    for path in cifs:
        array = load(path)
        seeds = C_HOT.mask(array, default="generate")
        print(f"--- {path.name} ---")
        print(f"  seed hotspots ({int(seeds.sum())}): {', '.join(describe(array, seeds))}")
        if not seeds.any():
            problems.append(f"{path.name}: no C_HOT seed to sample around")
            print("  FAIL: no C_HOT seed; the transform raises on this input")
            continue
        centroid = array.coord[seeds].astype(np.float64).mean(axis=0)
        spread = np.linalg.norm(array.coord[seeds].astype(np.float64) - centroid, axis=1)
        print(
            f"  seed offsets from their centroid: {', '.join(f'{d:.1f}' for d in spread)} A "
            f"({sampler.shape} patch radius {radius:.1f} A)"
        )

        np.random.seed(SEED)
        picked: collections.Counter = collections.Counter()
        size_hist: collections.Counter = collections.Counter()
        failures = 0
        seeds_dropped = 0
        for _ in range(args.n_draws):
            try:
                out = transform.forward({"atom_array": array, "example_id": path.stem})["atom_array"]
            except ValueError as exc:
                failures += 1
                if failures == 1:
                    print(f"  RAISED: {exc}")
                continue
            mask = C_HOT.mask(out, default="generate")
            if keeps_seeds and not bool((seeds & ~mask).sum() == 0):
                seeds_dropped += 1
            size_hist[int(mask.sum())] += 1
            picked.update(describe(out, mask & ~seeds))

        ok = args.n_draws - failures
        print(f"  draws: {ok}/{args.n_draws} succeeded, {failures} raised")
        if size_hist:
            print(f"  conditioned hotspots per design: {dict(sorted(size_hist.items()))}"
                  + (f"  (= {int(seeds.sum())} seeds + draw)" if keeps_seeds else ""))
        if keeps_seeds:
            kept = ok - seeds_dropped
            print(f"  seeds present in the final mask: {kept}/{ok} draws"
                  + ("" if seeds_dropped == 0 else "  <-- the wrapper is not doing its job"))
            if seeds_dropped:
                problems.append(f"{path.name}: {seeds_dropped}/{ok} draws lost a seed hotspot")
        print(f"  distinct SAMPLED atoms (seeds excluded): {len(picked)}")
        residues = sorted({lab.rsplit(" ", 1)[0] for lab in picked})
        print(f"  residues touched ({len(residues)}): {', '.join(residues) or '(none)'}")
        top = ", ".join(f"{lab} x{n}" for lab, n in picked.most_common(6))
        print(f"  most frequent sampled atoms: {top or '(none)'}")
        if failures:
            problems.append(f"{path.name}: {failures}/{args.n_draws} draws could not fit the requested count")
        if len(picked) < 3:
            problems.append(f"{path.name}: only {len(picked)} atom(s) are ever eligible; the patch is degenerate")
        print()

    print("=" * 70)
    if problems:
        print(f"{len(problems)} PROBLEM(S) -- do not run the patch arm with these settings:")
        for p in problems:
            print(f"  - {p}")
        print("Widen patch_diameter, lower min_distance/min_atom_sasa, or drop to num_hotspots=[2].")
        return 1
    print("Every draw fitted. This says nothing about design quality.")
    if keeps_seeds:
        print("The three seeds are conditioned on every design, plus the sampled atoms above.")
        print("hotspot_obedience is a FRACTION of ALL conditioned hotspots, so the campaign's")
        print("`min: 1.0` gate now demands every one of those 5-6 atoms within 4.5 A -- strictly")
        print("harder than the fixed arm's 3. Score the three seeds separately when aggregating.")
    else:
        print("The sampler excludes the seed RESIDUES, so A/204, A/206 and A/246 are never")
        print("themselves conditioned. Read the residue list above and decide whether that")
        print("patch is the epitope you want.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
