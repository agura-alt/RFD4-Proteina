"""Build annotated CIF inputs: CSS target (target_crop.pdb) + the four repo VHH frameworks.

Conditioning reproduces the scheme used by
``examples/antibody_design/structures/vhh/*.cif.zst`` and by the campaign template
``css_1g8pin_1G8_cdr123_L13-10-14_D166OD2_Q241NE2_Y246OH_generate.yaml``:

  target chain A  : absolute backbone coordinates (C_CRD, coordinate group 0)
                    full target sequence (C_SEQ)
                    three atom-level hotspots (C_HOT)
  binder chain H  : framework backbone in its own internal-coordinate group
                    (C_CRDGRP group 1 / C_INTCRD), dock NOT conditioned;
                    framework sequence known, CDR1/2/3 expandable sentinels free.
                    Carried over verbatim from the repo example file.

Hotspots are A/204 LEU CD2, A/206 LEU CD2, A/246 TYR CE1.
See the accompanying report for why CD2/CD2/CE1 and not CD1/CD1/CZ.

Run from the repository root.
"""

from __future__ import annotations

import contextlib
from pathlib import Path

import numpy as np

import rfproteina.conditions  # noqa: F401  -- registers condition classes before parse
from atomworks.io import parse
from atomworks.io.transforms.atom_array import (
    add_id_and_entity_annotations,
    add_iid_annotations,
)
from atomworks.io.utils.annotator import ensure_annotations
from atomworks.io.utils.atom_array_plus import (
    as_atom_array_plus,
    concatenate_atom_array_plus,
)
from atomworks.io.utils.io_utils import to_cif_file
from atomworks.ml.conditions import C_IDX, C_SEQ
from biotite.structure import BondType
from rfproteina.conditions import C_HOT
from rfproteina.conditions_reporting import (
    check_index_group_molecule_consistency,
    check_index_group_res_id_gaps,
    run_pipeline_validation,
    summarize_conditions,
    sync_index_groups,
)
from rfproteina.transforms.utils.conditioning_utils import (
    set_absolute_coordinate_conditioning,
)
from rfproteina.utils.parse_presets import ANNOTATED_CIF, DISTILLATION

ROOT = Path(__file__).resolve().parents[1]
TARGET_PDB = ROOT / "CSS" / "target_crop.pdb"
VHH_DIR = ROOT / "examples" / "antibody_design" / "structures" / "vhh"
OUT_DIR = ROOT / "CSS" / "inputs" / "cifs_L204CD2_L206CD2_Y246CE1"

FRAMEWORKS = [
    "VHH_caplacizumab",
    "VHH_ozoralizumab",
    "VHH_vobarilizumab",
    "h-NbBCII10",
]
# Chain H is atom-identical between the pmvk_3CH4 and mzb1_7AAH copies of each
# framework; only the (unconditioned) dock pose differs. pmvk is taken arbitrarily.
EXAMPLE_TARGET = "pmvk_3CH4"

TARGET_CHAIN = "A"
BINDER_CHAIN = "H"
BACKBONE = ("N", "CA", "C", "O")
HOTSPOTS = [(204, "LEU", "CD2"), (206, "LEU", "CD2"), (246, "TYR", "CE1")]
# Real disulfides, measured at 2.04 A; target_crop.pdb declares no SSBOND/CONECT
# and AtomWorks does not infer them, so they are added explicitly.
DISULFIDES = [(170, 195), (226, 264)]

PDB_CFG = DISTILLATION.replace(load_standard_annotations=False, return_atom_array_plus=True)

CHECKS = (
    run_pipeline_validation,
    check_index_group_res_id_gaps,
    check_index_group_molecule_consistency,
)


@contextlib.contextmanager
def _windows_randint_shim():
    """Work around a Windows-only bug in the design-validation pipeline.

    ``rfproteina/transforms/chirality_feature.py:171`` seeds with the legacy
    ``np.random.randint(0, 2**32)``. That legacy API defaults to dtype ``long``,
    which is int32 on Windows, so the call raises
    ``ValueError: high is out of bounds for int32`` and run_pipeline_validation
    fails before it validates anything. On Linux ``long`` is 64-bit and the call
    succeeds, which is why the repo (declared POSIX-only) has not hit this.

    The patched value is only used to seed a shuffle of chirality frames during
    validation; it does not touch the structure or any condition written to disk.
    """
    real = np.random.randint

    def patched(low, high=None, size=None, dtype=np.int64):
        return real(low, high, size, dtype)

    np.random.randint = patched
    try:
        yield
    finally:
        np.random.randint = real


def validate(aa, label: str) -> None:
    with _windows_randint_shim():
        violations = [msg for check in CHECKS for msg in check(aa)]
    if violations:
        raise ValueError(f"{label} validation failed:\n" + "\n".join(violations))
    print(f"    {label}: all three input-validation checks passed")


def load_target():
    aa = as_atom_array_plus(parse(str(TARGET_PDB), config=PDB_CFG)["assemblies"]["1"][0])
    ensure_annotations(aa, "is_protein", "is_ligand")

    # --- chemistry: add the two undeclared disulfides -------------------------
    for a, b in DISULFIDES:
        ia = np.flatnonzero((aa.res_id == a) & (aa.atom_name == "SG"))
        ib = np.flatnonzero((aa.res_id == b) & (aa.atom_name == "SG"))
        assert len(ia) == 1 and len(ib) == 1, f"SG lookup failed for CYS {a}/{b}"
        d = float(np.linalg.norm(aa.coord[ia[0]] - aa.coord[ib[0]]))
        assert d < 2.5, f"CYS{a}-CYS{b} SG-SG is {d:.2f} A; refusing to bond"
        aa.bonds.add_bond(int(ia[0]), int(ib[0]), BondType.SINGLE)
        print(f"    added disulfide CYS{a} SG - CYS{b} SG ({d:.2f} A)")

    target = aa.chain_id == TARGET_CHAIN
    assert target.all(), "target_crop.pdb is expected to hold chain A only"

    # --- conditioning ---------------------------------------------------------
    backbone = target & aa.is_protein & np.isin(aa.atom_name, BACKBONE)
    backbone &= np.isfinite(aa.coord).all(axis=1)
    assert backbone.sum() == 560, f"expected 560 backbone atoms, got {int(backbone.sum())}"
    set_absolute_coordinate_conditioning(aa, np.flatnonzero(backbone), True)

    C_SEQ.set_annotation_from_ground_truth(aa, target | aa.is_ligand)

    hot = np.zeros(aa.array_length(), dtype=bool)
    for rid, rname, aname in HOTSPOTS:
        m = target & (aa.res_id == rid) & (aa.atom_name == aname)
        assert m.sum() == 1, f"hotspot A/{rid}/{aname} matched {int(m.sum())} atoms"
        got = str(aa.res_name[m][0])
        assert got == rname, f"A/{rid} is {got}, expected {rname}"
        hot |= m
    C_HOT.set_annotation(aa, hot)
    assert hot.sum() == len(HOTSPOTS)

    C_IDX.set_annotation(aa, np.ones(aa.array_length(), dtype=bool))
    return aa


def load_framework(name: str):
    path = VHH_DIR / f"{EXAMPLE_TARGET}__{name}_cdr123.cif.zst"
    aa = as_atom_array_plus(parse(str(path), config=ANNOTATED_CIF)["assemblies"]["1"][0])
    binder = aa[aa.chain_id == BINDER_CHAIN]
    assert binder.array_length() > 0, f"no chain {BINDER_CHAIN} in {path.name}"
    return binder, path


def combine(target, binder):
    arrays = [target, binder]
    overlap = set(map(str, np.unique(target.chain_id))) & set(map(str, np.unique(binder.chain_id)))
    assert not overlap, f"chain id collision: {overlap}"

    # concatenate_atom_array_plus keeps only categories present in every input, so
    # materialize the union of condition/annotation columns on both first.
    columns = {
        name
        for arr in arrays
        for name in arr.get_annotation_categories()
        if name.startswith(("mask_", "condition_", "annotation_"))
    }
    for arr in arrays:
        ensure_annotations(arr, *columns)

    # The target comes from PDB and the framework from CIF, so a few plain structural
    # columns exist on only one side. Materialize them rather than dropping them:
    # auth_seq_id is the author numbering (identical to res_id for both inputs here),
    # and both alt-loc columns are '.' throughout since neither input has altlocs.
    for arr in arrays:
        ensure_annotations(arr, "is_protein", "is_ligand")
        cats = set(arr.get_annotation_categories())
        if "auth_seq_id" not in cats:
            arr.set_annotation("auth_seq_id", arr.res_id.astype(np.int64))
        for col in ("altloc_id", "label_alt_id"):
            if col not in cats:
                arr.set_annotation(col, np.full(arr.array_length(), ".", dtype="<U1"))

    combined = concatenate_atom_array_plus(arrays)
    combined = add_id_and_entity_annotations(combined, overwrite=True)
    combined = add_iid_annotations(combined, overwrite=True)
    sync_index_groups(combined, preserve_manual_edits=True)
    return combined


def compare_roundtrip(before, after, label: str) -> None:
    """Compare conditions atom-by-atom, keyed on atom identity, not on counts."""

    def key(aa):
        return list(
            zip(
                map(str, aa.chain_id),
                map(int, aa.res_id),
                map(str, aa.res_name),
                map(str, aa.atom_name),
            )
        )

    kb, ka = key(before), key(after)
    if kb != ka:
        raise ValueError(f"{label}: atom identity changed across the round trip")

    cols = sorted(
        {c for c in before.get_annotation_categories() if c.startswith(("mask_", "condition_", "annotation_"))}
        | {c for c in after.get_annotation_categories() if c.startswith(("mask_", "condition_", "annotation_"))}
    )
    for c in cols:
        if c not in before.get_annotation_categories() or c not in after.get_annotation_categories():
            raise ValueError(f"{label}: column {c} missing on one side of the round trip")
        vb, va = before.get_annotation(c), after.get_annotation(c)
        if np.issubdtype(vb.dtype, np.floating):
            ok = np.allclose(vb, va, equal_nan=True, atol=1e-4)
        else:
            ok = bool(np.array_equal(vb, va))
        if not ok:
            raise ValueError(f"{label}: condition column {c} changed across the round trip")

    # chemistry
    bb = before.bonds.as_array() if before.bonds is not None else np.empty((0, 3), int)
    ba = after.bonds.as_array() if after.bonds is not None else np.empty((0, 3), int)
    sb = {(min(i, j), max(i, j), t) for i, j, t in bb.tolist()}
    sa = {(min(i, j), max(i, j), t) for i, j, t in ba.tolist()}
    if sb != sa:
        raise ValueError(f"{label}: bond set changed across the round trip ({len(sb)} -> {len(sa)})")
    if not np.array_equal(before.charge, after.charge):
        raise ValueError(f"{label}: formal charges changed across the round trip")
    print(f"    round trip: {len(cols)} condition columns, {len(sb)} bonds, charges all identical")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading target (CSS/target_crop.pdb, preset DISTILLATION, no missing-atom repair)")
    target = load_target()
    n_charged = int((target.charge != 0).sum())
    print(f"    {target.array_length()} atoms, {len(np.unique(target.res_id))} residues, "
          f"{n_charged} formally charged atoms, net {int(target.charge.sum()):+d}")

    for fw in FRAMEWORKS:
        print(f"\n=== {fw} ===")
        binder, src = load_framework(fw)
        print(f"    framework from {src.name}: {binder.array_length()} atoms, "
              f"{len(np.unique(binder.res_id))} residues")

        combined = combine(target.copy(), binder)
        print(f"    combined: {combined.array_length()} atoms, chains {list(np.unique(combined.chain_id))}")
        validate(combined, "pre-save")

        out = OUT_DIR / f"css_{fw}_cdr123_L204CD2_L206CD2_Y246CE1.cif"
        to_cif_file(
            combined,
            str(out),
            save_standard_annotations=True,
            chain_disambiguation="chain_iid",
        )
        reloaded = as_atom_array_plus(
            parse(str(out), config=ANNOTATED_CIF)["assemblies"]["1"][0]
        )
        validate(reloaded, "round-trip")
        compare_roundtrip(combined, reloaded, fw)

        print("    --- conditions on the saved file ---")
        for line in summarize_conditions(reloaded).splitlines():
            if line.strip():
                print("      " + line)
        print(f"    wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
