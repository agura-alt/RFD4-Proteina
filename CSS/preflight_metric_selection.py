"""Offline preflight: which metrics will actually run, given metrics_tags?

Composes the campaign config with Hydra and applies the same tag rule the runtime
uses (rfproteina/metrics/metric.py:415) to every configured metric, without a GPU,
a checkpoint, or an AF3 call. Answers "did adding the `antibody` tag really turn
anything on" before you spend a smoke run finding out.

Usage, from the repository root:
    python CSS/preflight_metric_selection.py
"""

from __future__ import annotations

import os
from pathlib import Path

from hydra import compose, initialize_config_dir
from omegaconf import OmegaConf

ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "CSS" / "configs"
CONFIG_NAME = "css_cdr123_L204CD2_L206CD2_Y246CE1_generate"

# Tags the runtime adds per refold; included so refolding metrics resolve as they
# would once both oracles have run.
RUNTIME_REFOLD_TAGS = {
    "is_refold:oracle_ppi",
    "exists_refold:oracle_ppi",
    "is_refold:oracle_apo.apo",
    "exists_refold:oracle_apo.apo",
}

# Metrics the campaign is counting on. Flagged loudly if they resolve to "off".
EXPECTED_ON = [
    "hotspot_obedience",
    "interchain_interactions",
    "antibody_contacts",
    "design_motif_allatom_rmsd",
    "sequence_adherence",
]


def matches(entry, tags: set[str]) -> bool:
    """Mirror Metric.matches_tags: metric declares, tags are intersected against it."""
    req_all = set(entry.get("required_tags_all") or [])
    req_any = set(entry.get("required_tags_any") or [])
    prohibited = set(entry.get("prohibited_tags") or [])
    return (
        (not req_all or tags.issuperset(req_all))
        and (not req_any or bool(tags & req_any))
        and not (tags & prohibited)
    )


def main() -> int:
    os.environ.setdefault("PROJECT_ROOT", str(ROOT))
    os.environ.setdefault("DATA_PATH", "/unused/for/preflight")

    with initialize_config_dir(config_dir=str(CONFIG_DIR), version_base=None):
        cfg = compose(config_name=CONFIG_NAME)

    tags = set(OmegaConf.to_container(cfg.metrics_tags))
    print("metrics_tags:")
    for t in sorted(tags):
        print(f"    {t}")

    # --- recycling ---------------------------------------------------------
    gen = cfg.generation
    schedule = OmegaConf.to_container(gen.n_recycle_schedule) or []
    steps = sum(n for n, _ in schedule)
    passes = sum(n * k for n, k in schedule)
    print(f"\nrecycling: {steps} sampling steps, {passes} total trunk passes")
    print(f"    schedule (steps, recycles): {[tuple(x) for x in schedule]}")
    print(f"    generation.args.nsteps = {gen.args.nsteps}")
    if passes == 180 and steps == 60:
        print("    ok -- matches PPI medium/front_loaded from the docs table")
    else:
        print("    WARNING -- does not match medium/front_loaded (60 steps / 180 passes)")
    if not schedule:
        print("    FAIL -- empty schedule means zero recycling; the designs will not dock")

    # --- metric selection --------------------------------------------------
    full = tags | RUNTIME_REFOLD_TAGS
    groups = ["inline_metrics", "refolding_metrics", "refolding_oracles", "pass_criteria"]
    selected: dict[str, list[str]] = {}
    skipped: dict[str, list[str]] = {}

    for group in groups:
        node = cfg.metrics.get(group)
        selected[group], skipped[group] = [], []
        if not node:
            continue
        for name, entry in node.items():
            if not isinstance(entry, (dict, type(cfg))) or entry is None:
                continue
            e = OmegaConf.to_container(entry, resolve=False)
            if not isinstance(e, dict):
                continue
            (selected if matches(e, full) else skipped)[group].append(name)

    for group in groups:
        print(f"\n--- {group} ---")
        print(f"  ON  ({len(selected[group])}): {', '.join(sorted(selected[group])) or '(none)'}")
        print(f"  off ({len(skipped[group])}): {', '.join(sorted(skipped[group])) or '(none)'}")

    all_on = {n for g in groups for n in selected[g]}
    print("\n--- metrics this campaign depends on ---")
    problems = []
    for name in EXPECTED_ON:
        if name in all_on:
            print(f"  ok      {name}")
        else:
            print(f"  OFF     {name}")
            problems.append(name)

    print("\n" + "=" * 60)
    if problems:
        print("These expected metrics resolve to OFF: " + ", ".join(problems))
        return 1
    print("Metric selection looks right. This is a static tag check only --")
    print("a metric can still emit no column at runtime if its masks select 0 atoms.")
    print("CSS/audit_smoke_metrics.py is what checks that, after the smoke run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
