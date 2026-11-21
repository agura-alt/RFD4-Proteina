"""Build the "high constraint" arm: CSS target (target_crop.pdb) + the four repo VHH frameworks.

Template dir: CSS/inputs/cifs_L204CD2_L206CD2_Y246CE1/
Everything is reproduced from that build (CSS/build_annotated_cifs.py) except C_HOT.
Each output is diffed against its own framework's template; only
``condition_hotspot_1_atom`` is permitted to differ.

Hotspots (six): A/165 GLU OE2, A/166 ASP OD2, A/204 LEU CD2, A/206 LEU CD2,
                A/241 GLN NE2, A/246 TYR CE1.

Atom provenance -- "the same atom that was used in the previous run":
  204 CD2, 206 CD2, 246 CE1 -- the run named as the template (cifs_L204CD2_L206CD2_Y246CE1).
  166 OD2, 241 NE2          -- the 1G8 positive control
                               (css_1g8pin_1G8_cdr123_L13-10-14_D166OD2_Q241NE2_Y246OH_generate.yaml).
  165 OE2                   -- the "spanning the groove" arm
                               (cifs_E165OE2_L204CD2_L206CD2_Q241NE2), user-confirmed 2026-11-19.
See the accompanying REPORT.md for the TYR 246 CE1-vs-OH conflict between those two runs.

Run from the repository root.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "CSS"))

import rfproteina.conditions  # noqa: F401  -- registers condition classes before parse
from atomworks.io import parse
from atomworks.io.utils.atom_array_plus import as_atom_array_plus
from atomworks.io.utils.io_utils import to_cif_file
from rfproteina.conditions import C_HOT
from rfproteina.conditions_reporting import summarize_conditions
from rfproteina.utils.parse_presets import ANNOTATED_CIF

# Reuse the verified machinery from the previous builds verbatim.
import build_annotated_cifs as prev

FRAMEWORKS = [
    "VHH_caplacizumab",
    "VHH_ozoralizumab",
    "VHH_vobarilizumab",
    "h-NbBCII10",
]
TAG = "E165OE2_D166OD2_L204CD2_L206CD2_Q241NE2_Y246CE1"
TEMPLATE_DIR = ROOT / "CSS" / "inputs" / "cifs_L204CD2_L206CD2_Y246CE1"
OUT_DIR = ROOT / "CSS" / "inputs" / f"cifs_{TAG}"

HOTSPOTS = [
    (165, "GLU", "OE2"),
    (166, "ASP", "OD2"),
    (204, "LEU", "CD2"),
    (206, "LEU", "CD2"),
    (241, "GLN", "NE2"),
    (246, "TYR", "CE1"),
]


def key(aa):
    return list(zip(map(str, aa.chain_id), map(int, aa.res_id),
                    map(str, aa.res_name), map(str, aa.atom_name)))


def condition_cols(aa):
    return {c for c in aa.get_annotation_categories()
            if c.startswith(("mask_", "condition_", "annotation_"))}


def diff_against_template(reloaded, template_path: Path) -> None:
    tmpl = as_atom_array_plus(parse(str(template_path), config=ANNOTATED_CIF)["assemblies"]["1"][0])
    if tmpl.array_length() != reloaded.array_length():
        raise ValueError(f"atom count differs: {tmpl.array_length()} vs {reloaded.array_length()}")
    if key(tmpl) != key(reloaded):
        raise ValueError("atom identity differs from the template")
    if not np.allclose(tmpl.coord, reloaded.coord, equal_nan=True, atol=1e-3):
        raise ValueError("coordinates differ from the template")
    if not np.array_equal(tmpl.charge, reloaded.charge):
        raise ValueError("formal charges differ from the template")
    sb = {(min(i, j), max(i, j), t) for i, j, t in tmpl.bonds.as_array().tolist()}
    sa = {(min(i, j), max(i, j), t) for i, j, t in reloaded.bonds.as_array().tolist()}
    if sb != sa:
        raise ValueError(f"bond set differs from the template ({len(sb)} vs {len(sa)})")

    cols = sorted(condition_cols(tmpl) | condition_cols(reloaded))
    changed = []
    for c in cols:
        vb, va = tmpl.get_annotation(c), reloaded.get_annotation(c)
        same = (np.allclose(vb, va, equal_nan=True, atol=1e-4)
                if np.issubdtype(vb.dtype, np.floating) else bool(np.array_equal(vb, va)))
        if not same:
            changed.append(c)
    expected = [c for c in cols if "hot" in c.lower()]
    if sorted(changed) != sorted(expected):
        raise ValueError(f"expected only {expected} to change, got {changed}")

    th = C_HOT.mask(tmpl, default="generate")
    print(f"    vs template: identity/coords/{len(sa)} bonds/charges identical; "
          f"{len(cols) - len(changed)}/{len(cols)} condition columns identical; "
          f"only {changed[0]} differs")
    print("      template C_HOT: " + ", ".join(
        f"A/{int(r)}/{n}/{a}" for r, n, a in zip(tmpl.res_id[th], tmpl.res_name[th], tmpl.atom_name[th])))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading target (CSS/target_crop.pdb, DISTILLATION, no missing-atom repair)")
    prev.HOTSPOTS = HOTSPOTS  # load_target() reads this module-level list
    target = prev.load_target()
    print(f"    {target.array_length()} atoms, {len(np.unique(target.res_id))} residues, "
          f"net charge {int(target.charge.sum()):+d}")
    hot = C_HOT.mask(target, default="generate")
    print(f"    C_HOT: {int(hot.sum())} atoms -> "
          + ", ".join(f"A/{r}/{n}/{a}" for r, n, a in HOTSPOTS))

    for fw in FRAMEWORKS:
        print(f"\n=== {fw} ===")
        template = TEMPLATE_DIR / f"css_{fw}_cdr123_L204CD2_L206CD2_Y246CE1.cif"
        assert template.exists(), f"template not found: {template}"

        binder, src = prev.load_framework(fw)
        print(f"    framework chain H from {src.name}: {binder.array_length()} atoms")

        combined = prev.combine(target.copy(), binder)
        print(f"    combined: {combined.array_length()} atoms, "
              f"chains {list(np.unique(combined.chain_id))}")
        prev.validate(combined, "pre-save")

        out = OUT_DIR / f"css_{fw}_cdr123_{TAG}.cif"
        to_cif_file(combined, str(out), save_standard_annotations=True,
                    chain_disambiguation="chain_iid")
        reloaded = as_atom_array_plus(parse(str(out), config=ANNOTATED_CIF)["assemblies"]["1"][0])
        prev.validate(reloaded, "round-trip")
        prev.compare_roundtrip(combined, reloaded, fw)
        diff_against_template(reloaded, template)

        rh = C_HOT.mask(reloaded, default="generate")
        print("      new      C_HOT: " + ", ".join(
            f"A/{int(r)}/{n}/{a}" for r, n, a in
            zip(reloaded.res_id[rh], reloaded.res_name[rh], reloaded.atom_name[rh])))
        print(f"    wrote {out.relative_to(ROOT)}")

    print("\n--- conditions on the last saved file ---")
    for line in summarize_conditions(reloaded).splitlines():
        if line.strip():
            print("  " + line)


if __name__ == "__main__":
    main()
