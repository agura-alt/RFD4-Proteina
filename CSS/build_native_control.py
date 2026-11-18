"""Build the native 1G8 complex as a positive control for `oracle_ppi`.

Why this exists
---------------
The campaign's designs score `oracle_ppi.min_ipae_coi` of 15.8-24.5 against a gate of 1.2
(smoke run 2026-11-18, 12 designs). That number is uninterpretable on its own: nobody has
measured what this oracle reports for a binder that is KNOWN to work on this antigen, so
there is no way to tell a failing design from an unreachable gate.

`1G8-PTRAMP-CSS.pdb` is a solved structure of a real nanobody bound to this antigen. Scoring
that native complex through the same oracle gives the ceiling. Interpretation:

  native ~<=1    the gate is reachable here; the designs genuinely do not dock.
  native 10-25   a crystallographically verified binder fails the gate, so the gate is not
                 measuring design quality on this target and every campaign number that
                 rests on it needs reinterpreting.

This script only BUILDS the input and reports whether the control is meaningful. It runs no
oracle; the `rfd4-evaluate` command is printed at the end.

Caveats this script checks for you
----------------------------------
* The PDB holds two copies of the complex (A/B/C and D/E/F). Chains C and F are the
  nanobodies -- they carry the canonical VHH disulfide CYS22-CYS92 (SSBOND 7 and 14).
* The campaign target is chain A cropped to residues 151-290. 1G8's own epitope was chosen
  separately from this campaign's (L204 CD2 / L206 CD2 / Y246 CE1), so the native nanobody
  could in principle bind outside the crop. Measured 2026-11-18, it does not: chain C buries
  21 of the 140 cropped residues, and reaches all three conditioned hotspot atoms --
      A/204 CD2  3.02 A      A/206 CD2  3.55 A      A/246 CE1  3.61 A
  all inside the 4.5 A `max_hotspot_distance` that `hotspot_adherence` itself uses. The
  native binder would pass this campaign's hotspot gate. Its interface also spans the 1G8
  arm's D166/Q241/Y246 patch, so both epitopes are faces of one real paratope. That makes
  this a direct positive control at this campaign's own epitope, not merely a calibration of
  the metric on this antigen. The script re-measures all of it rather than trusting this note.

Usage, from the repository root:
    python CSS/build_native_control.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE_PDB = ROOT / "CSS" / "1G8-PTRAMP-CSS.pdb"
OUT_DIR = ROOT / "CSS" / "native_control"
OUT_CIF = OUT_DIR / "1g8_native_complex.cif"

TARGET_CHAIN = "A"
TARGET_FIRST, TARGET_LAST = 151, 290  # matches CSS/target_crop.pdb
NANOBODY_CHAIN = "C"  # paired copy of the complex; F is the second
BINDER_CHAIN_OUT = "H"  # campaign convention; the config's masks address the binder as H

# The campaign's three conditioned hotspot atoms, in source numbering.
HOTSPOTS = [("LEU", 204, "CD2"), ("LEU", 206, "CD2"), ("TYR", 246, "CE1")]

CONTACT_CUTOFF = 5.0  # A, heavy-atom; generous, this is a "does it touch at all" check


def _min_distance(a: struc.AtomArray, b: struc.AtomArray) -> float:
    """Smallest heavy-atom distance between two AtomArrays."""
    if a.array_length() == 0 or b.array_length() == 0:
        return float("inf")
    diff = a.coord[:, None, :] - b.coord[None, :, :]
    return float(np.sqrt((diff**2).sum(-1)).min())


def main() -> int:
    if not SOURCE_PDB.is_file():
        print(f"FAIL: {SOURCE_PDB} not found")
        return 1

    aa = pdb.get_structure(pdb.PDBFile.read(str(SOURCE_PDB)), model=1, altloc="first")
    # Drop waters, ions and other heteroatoms; the oracle folds polymers.
    aa = aa[struc.filter_amino_acids(aa) & ~aa.hetero]

    present = set(np.unique(aa.chain_id))
    for chain in (TARGET_CHAIN, NANOBODY_CHAIN):
        if chain not in present:
            print(f"FAIL: chain {chain!r} absent from {SOURCE_PDB.name}; found {sorted(present)}")
            return 1

    target = aa[
        (aa.chain_id == TARGET_CHAIN) & (aa.res_id >= TARGET_FIRST) & (aa.res_id <= TARGET_LAST)
    ].copy()
    nanobody = aa[aa.chain_id == NANOBODY_CHAIN].copy()

    if target.array_length() == 0 or nanobody.array_length() == 0:
        print("FAIL: target crop or nanobody selection is empty")
        return 1

    print("--- chains ---")
    t_res = np.unique(target.res_id)
    n_res = np.unique(nanobody.res_id)
    print(f"  target   chain {TARGET_CHAIN} {t_res.min()}-{t_res.max()}  "
          f"{len(t_res)} residues, {target.array_length()} atoms")
    print(f"  nanobody chain {NANOBODY_CHAIN} {n_res.min()}-{n_res.max()}  "
          f"{len(n_res)} residues, {nanobody.array_length()} atoms")

    # Sanity: the canonical VHH disulfide should be present and ~2.05 A.
    sg = nanobody[(nanobody.atom_name == "SG") & (nanobody.res_name == "CYS")]
    if sg.array_length() >= 2:
        pairs = [
            (int(sg.res_id[i]), int(sg.res_id[j]), float(np.linalg.norm(sg.coord[i] - sg.coord[j])))
            for i in range(sg.array_length())
            for j in range(i + 1, sg.array_length())
        ]
        bonded = [p for p in pairs if p[2] < 2.5]
        print(f"  nanobody disulfides: {[(a, b, round(d, 2)) for a, b, d in bonded] or 'none under 2.5 A'}")

    # --- is this control meaningful? ------------------------------------
    print(f"\n--- does the native nanobody contact the crop? (heavy atom, {CONTACT_CUTOFF} A) ---")
    contacting = []
    for res_id in t_res:
        res_atoms = target[target.res_id == res_id]
        d = _min_distance(res_atoms, nanobody)
        if d <= CONTACT_CUTOFF:
            contacting.append((int(res_id), str(res_atoms.res_name[0]), round(d, 2)))
    print(f"  {len(contacting)} of {len(t_res)} cropped target residues within {CONTACT_CUTOFF} A")
    if contacting:
        shown = ", ".join(f"{name}{rid}" for rid, name, _ in contacting[:20])
        print(f"  {shown}{' ...' if len(contacting) > 20 else ''}")

    print("\n--- campaign hotspot atoms vs the native nanobody ---")
    for res_name, res_id, atom_name in HOTSPOTS:
        sel = target[
            (target.res_id == res_id) & (target.res_name == res_name) & (target.atom_name == atom_name)
        ]
        if sel.array_length() == 0:
            print(f"  {res_name}{res_id} {atom_name:4s} NOT FOUND in the crop")
            continue
        d = _min_distance(sel, nanobody)
        verdict = "contact" if d <= CONTACT_CUTOFF else "no contact"
        print(f"  {res_name}{res_id} {atom_name:4s} nearest nanobody atom {d:6.2f} A  ({verdict})")

    if not contacting:
        print(
            "\nWARNING: the native nanobody does not touch the cropped target at all.\n"
            "  It binds elsewhere on chain A, or binds chain B. The refold would then score a\n"
            "  non-complex and tell you nothing. Either widen the crop to 1G8's real epitope or\n"
            "  use the full chain A before trusting the number."
        )

    # --- write the control input ----------------------------------------
    nanobody.chain_id = np.full(nanobody.array_length(), BINDER_CHAIN_OUT)
    complex_aa = target + nanobody

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cif = pdbx.CIFFile()
    pdbx.set_structure(cif, complex_aa)
    cif.write(str(OUT_CIF))

    print(f"\nwrote {OUT_CIF}")
    print(f"  chains {sorted(set(np.unique(complex_aa.chain_id)))}, {complex_aa.array_length()} atoms")
    print(
        "\nScore it with the campaign's own tags:\n"
        f"  rfd4-evaluate inputs={OUT_DIR} metrics_tags=[ppi_refolding,antibody] \\\n"
        f"      out_dir=./runs/native_control verbose=true\n"
        "\nThen read `oracle_ppi.min_ipae_A_1_H_1` from\n"
        "  ./runs/native_control/evaluation_metrics/per_example_metrics.csv\n"
        "The `*_coi` columns will be empty: chain-of-interest is picked as the least\n"
        "sequence-conditioned chain, and on an unannotated crystal structure both chains tie,\n"
        "so it resolves as ambiguous. The pairwise column is computed regardless and is the\n"
        "same quantity (identical to min_ipae_coi on every design of the 2026-11-18 smoke run).\n"
        "\nCompare against the designs' 15.8-24.5."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
