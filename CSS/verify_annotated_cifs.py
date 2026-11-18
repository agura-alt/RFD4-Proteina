"""Verify the built annotated CIFs against their two sources, atom by atom.

Checks, per output file:
  * chain A is target_crop.pdb exactly (atom identity, coordinates, formal charges)
  * chain A carries the two disulfides that the PDB left undeclared
  * chain H is the repo example's framework chain verbatim (identity, charges,
    internal-coordinate conditioning, sequence conditioning, CDR sentinels)
  * conditioning is the intended one and nothing else is conditioned

Run from the repository root.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import rfproteina.conditions  # noqa: F401
from atomworks.io import parse
from atomworks.io.utils.atom_array_plus import as_atom_array_plus
from rfproteina.utils.parse_presets import ANNOTATED_CIF, DISTILLATION

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "CSS" / "inputs" / "cifs_L204CD2_L206CD2_Y246CE1"
VHH_DIR = ROOT / "examples" / "antibody_design" / "structures" / "vhh"
FRAMEWORKS = ["VHH_caplacizumab", "VHH_ozoralizumab", "VHH_vobarilizumab", "h-NbBCII10"]
HOTSPOTS = {(204, "CD2"), (206, "CD2"), (246, "CE1")}
PDB_CFG = DISTILLATION.replace(load_standard_annotations=False, return_atom_array_plus=True)

failures: list[str] = []


def check(cond: bool, msg: str) -> None:
    if cond:
        print(f"      ok   {msg}")
    else:
        print(f"      FAIL {msg}")
        failures.append(msg)


def ident(aa):
    return list(zip(map(int, aa.res_id), map(str, aa.res_name), map(str, aa.atom_name)))


def disulfides(aa):
    if aa.bonds is None:
        return set()
    sg = set(np.flatnonzero(aa.atom_name == "SG").tolist())
    return {
        tuple(sorted((int(aa.res_id[i]), int(aa.res_id[j]))))
        for i, j, _ in aa.bonds.as_array().tolist()
        if i in sg and j in sg
    }


print("Loading sources")
target_src = as_atom_array_plus(parse(str(ROOT / "CSS/target_crop.pdb"), config=PDB_CFG)["assemblies"]["1"][0])
print(f"  target_crop.pdb: {target_src.array_length()} atoms")

for fw in FRAMEWORKS:
    print(f"\n=== {fw} ===")
    out = OUT_DIR / f"css_{fw}_cdr123_L204CD2_L206CD2_Y246CE1.cif"
    built = as_atom_array_plus(parse(str(out), config=ANNOTATED_CIF)["assemblies"]["1"][0])
    ex = as_atom_array_plus(
        parse(str(VHH_DIR / f"pmvk_3CH4__{fw}_cdr123.cif.zst"), config=ANNOTATED_CIF)["assemblies"]["1"][0]
    )
    fw_src = ex[ex.chain_id == "H"]

    A = built[built.chain_id == "A"]
    H = built[built.chain_id == "H"]

    print("    -- target chain A vs target_crop.pdb --")
    check(ident(A) == ident(target_src), "atom identity identical (1134 atoms, res 151-290)")
    check(np.allclose(A.coord, target_src.coord, atol=1e-3), "coordinates identical")
    check(np.array_equal(A.charge, target_src.charge), f"formal charges identical ({int((A.charge != 0).sum())} charged, net {int(A.charge.sum()):+d})")
    check(disulfides(A) == {(170, 195), (226, 264)}, f"disulfides present: {sorted(disulfides(A))}")

    print("    -- binder chain H vs repo example --")
    check(ident(H) == ident(fw_src), f"atom identity identical ({H.array_length()} atoms)")
    check(np.array_equal(H.charge, fw_src.charge), f"formal charges identical ({int((H.charge != 0).sum())} charged)")
    for col in ("condition_InternalCoordinate_1_atom", "condition_CoordinateGroup_1_atom",
                "condition_sequence_1_residue", "condition_cath-label_1_residue",
                "annotation_expsegmin_1_residue", "annotation_expsegmax_1_residue"):
        a, b = H.get_annotation(col), fw_src.get_annotation(col)
        same = np.allclose(a, b, equal_nan=True) if np.issubdtype(a.dtype, np.floating) else np.array_equal(a, b)
        check(bool(same), f"{col} carried over verbatim")
    segs = sorted({(int(lo), int(hi)) for lo, hi in zip(
        H.get_annotation("annotation_expsegmin_1_residue"),
        H.get_annotation("annotation_expsegmax_1_residue")) if hi > 0})
    check(segs == [(9, 10), (12, 21), (13, 13)], f"CDR expandable segments intact: {segs}")

    print("    -- conditioning --")
    grp = built.get_annotation("condition_CoordinateGroup_1_atom")
    crd = built.get_annotation("condition_coordinate_1_atom")
    icd = built.get_annotation("condition_InternalCoordinate_1_atom")
    hot = built.get_annotation("condition_hotspot_1_atom")
    bbA = (built.chain_id == "A") & np.isin(built.atom_name, ["N", "CA", "C", "O"])

    check(int((grp == 0).sum()) == 560 and bool(np.array_equal(grp == 0, bbA)),
          "coordinate group 0 == target backbone exactly (560 atoms, absolute)")
    check(bool(np.array_equal(np.isfinite(crd).all(axis=1), grp == 0)),
          "absolute coordinate targets exist exactly on group 0")
    check(int((grp == 1).sum()) == 332 and bool(((built.chain_id == "H") | (grp != 1)).all()),
          "coordinate group 1 == framework backbone only (332 atoms, internal)")
    check(not np.isfinite(crd[built.chain_id == "H"]).any(),
          "framework dock NOT conditioned (absolute coords all NaN on chain H)")
    check(bool(np.array_equal(np.isfinite(icd).all(axis=1), grp == 1)),
          "internal coordinate targets exist exactly on group 1")
    hotset = {(int(built.res_id[i]), str(built.atom_name[i])) for i in np.flatnonzero(hot)}
    check(hotset == HOTSPOTS and bool((built.chain_id[hot] == "A").all()),
          f"hotspots are exactly the intended target atoms: {sorted(hotset)}")

    seq = built.get_annotation("condition_sequence_1_residue")
    unknown_A = (built.chain_id == "A") & (seq == "<M>")
    check(not unknown_A.any(), "full target sequence is conditioned (no masked target residues)")
    unknown_H = (built.chain_id == "H") & (seq == "<M>")
    check(int(unknown_H.sum()) == 3, "framework sequence known; the 3 CDR sentinels left free")

print("\n" + "=" * 60)
if failures:
    print(f"{len(failures)} CHECK(S) FAILED")
    for f in failures:
        print("  -", f)
    raise SystemExit(1)
print("All verification checks passed for all four files.")
