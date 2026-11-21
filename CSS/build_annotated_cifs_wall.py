"""Build the "wall only hydrophilic" arm: CSS target (target_crop.pdb) + the four repo VHH frameworks.

"wall only hydrophilic" is a colloquial arm name only; it carries no semantics for the
build. All four hotspot residues are non-hydrophobic in CSS_epitope_residue_report.csv,
and the two leucines of the previous arms are dropped, so the name is at least accurate.

Template dir: CSS/inputs/cifs_L204CD2_L206CD2_Y246CE1/
Everything is reproduced from that build (CSS/build_annotated_cifs.py) except C_HOT.
Each output is diffed against its own framework's template; only
``condition_hotspot_1_atom`` is permitted to differ.

Hotspots: A/165 GLU OE2, A/166 ASP OD2, A/239 HIS CE1, A/241 GLN NE2.
  165 OE2 -- the atom used in the previous (E165OE2_L204CD2_L206CD2_Q241NE2) run.
  166 OD2 -- the atom used in the 1G8 positive-control run (..._D166OD2_...).
  241 NE2 -- the atom used in the 1G8 positive control and in the E165OE2 run.
  239 CE1 -- NO previous run used residue 239. CE1 is the most solvent-exposed atom
             of the imidazole (26.93 A^2). Chosen by the user 2026-11-21 over NE2
             (18.26 A^2, the most exposed *nitrogen*), which the first draft of this
             arm used; the switch takes raw exposure over H-bond capability, so the
             "hydrophilic" name no longer describes this contact -- CE1 is a ring
             carbon. ND1 is buried (2.52 A^2) and would have been the HIS analogue of
             the LEU CD1 trap resolved in the first build. See REPORT.md.

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

import build_annotated_cifs as prev
# Reuse the diff/compare machinery of the groove build verbatim.
import build_annotated_cifs_groove as groove

FRAMEWORKS = groove.FRAMEWORKS
TEMPLATE_DIR = ROOT / "CSS" / "inputs" / "cifs_L204CD2_L206CD2_Y246CE1"
OUT_DIR = ROOT / "CSS" / "inputs" / "cifs_E165OE2_D166OD2_H239CE1_Q241NE2"

TAG = "E165OE2_D166OD2_H239CE1_Q241NE2"

HOTSPOTS = [
    (165, "GLU", "OE2"),
    (166, "ASP", "OD2"),
    (239, "HIS", "CE1"),
    (241, "GLN", "NE2"),
]


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
        groove.diff_against_template(reloaded, template)

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
