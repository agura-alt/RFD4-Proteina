"""Build the native 1G8 complex as an ANNOTATED CIF, so the campaign's own metrics run on it.

Why a second builder
--------------------
`build_native_control.py` writes a plain two-chain structure with biotite. That is enough for
a bare AF3 fold, and it produced the headline number (native `min_ipae` 16.56, inside the
designs' 15.8-24.5), but it carries no conditioning, so:

  * chain-of-interest cannot be resolved -- `annotate_chain_of_interest` raises, which takes
    the whole refolding pipeline down with it and leaves no metrics at all;
  * `oracle_ppi`'s template query finds no `C_CRD` and silently templates nothing, so the
    control is scored with less information than the designs were;
  * `oracle_apo` has no designed chain to isolate, so `apo_holo_consistency` cannot run.

This builder conditions the complex the way a design input is conditioned, so the native
binder goes through exactly the machinery the designs did. It reuses `load_target()` from
`build_annotated_cifs.py` verbatim rather than reimplementing it, so the target half is
identical by construction -- same crop, same two disulfides, same hotspots, same C_CRD/C_SEQ.

What the binder gets
--------------------
  C_SEQ on framework positions only, absent across CDR1/2/3.

That mirrors a design, where the framework sequence was given and the loops were generated,
and it is what `refolded_apo_holo_cdr_rmsd` keys on: it aligns on `sa_mask('C_SEQ')` and
measures the complement. It also makes the binder the least sequence-conditioned chain, so
chain-of-interest resolves to H and `ApoOnlyTransform` isolates the right chain.

Deliberately NOT set: internal-coordinate conditioning on the binder. A design input carries
it to tell the generator to keep framework geometry, but nothing in the evaluation path reads
it, and fabricating it here would add a large guess for no gain. This file is an input to
`rfd4-evaluate`, not to generation; the validation checks below are reported rather than
enforced for that reason.

CDR boundaries -- CHECK THESE
-----------------------------
Chain C is antibody-numbered with insertion codes (52A, 82A-C, 100A-D), so CDRs are taken by
number rather than guessed from geometry. The ranges below are the usual Kabat ones and the
script prints the residues and lengths it selected. The 1G8 config states this nanobody's
native loop lengths are 13 / 10 / 14; if the printed lengths disagree, adjust `CDR_RANGES`
rather than trusting the defaults. The whole-binder metric `refolded_apo_holo_rmsd` does not
depend on this -- only the CDR-specific one does.

Usage, from the repository root:
    python CSS/build_native_control_annotated.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# Same directory; importing executes constants and helpers but not main().
sys.path.insert(0, str(Path(__file__).resolve().parent))

import rfproteina.conditions  # noqa: F401  -- registers condition classes before parse
from atomworks.io import parse
from atomworks.io.utils.annotator import ensure_annotations
from atomworks.io.utils.atom_array_plus import as_atom_array_plus
from atomworks.io.utils.io_utils import to_cif_file
from atomworks.ml.conditions import C_IDX, C_SEQ
from biotite.structure import BondType, get_residue_starts
from rfproteina.conditions_reporting import summarize_conditions
from rfproteina.utils.parse_presets import ANNOTATED_CIF

from build_annotated_cifs import (  # noqa: E402
    PDB_CFG,
    ROOT,
    combine,
    compare_roundtrip,
    load_target,
    validate,
)

SOURCE_PDB = ROOT / "CSS" / "1G8-PTRAMP-CSS.pdb"
OUT_DIR = ROOT / "CSS" / "inputs" / "native_control"
OUT_CIF = OUT_DIR / "css_native_1G8_control.cif"

NANOBODY_CHAIN = "C"  # the second copy, F, is atom-identical
BINDER_CHAIN = "H"

# Kabat CDR ranges in the source antibody numbering. Insertion codes inside a range are
# included automatically, which is how 52A lands in CDR2 and 100A-D in CDR3.
CDR_RANGES = {"CDR1": (31, 35), "CDR2": (50, 65), "CDR3": (95, 102)}
EXPECTED_CDR_LENGTHS = (13, 10, 14)  # per the 1G8 config; advisory only

# 1G8's canonical VHH disulfide. The target builder adds its two explicitly because
# AtomWorks infers none; the binder needs the same treatment.
BINDER_DISULFIDE = (22, 92)


def load_native_binder():
    """Chain C as a conditioned binder: CDR masks taken in antibody numbering, then renumbered."""
    aa = as_atom_array_plus(parse(str(SOURCE_PDB), config=PDB_CFG)["assemblies"]["1"][0])
    ensure_annotations(aa, "is_protein", "is_ligand")
    binder = aa[aa.chain_id == NANOBODY_CHAIN].copy()
    if binder.array_length() == 0:
        raise ValueError(f"no chain {NANOBODY_CHAIN!r} in {SOURCE_PDB.name}")

    # --- CDR masks, while antibody numbering is still intact ------------------
    ins = (
        binder.get_annotation("ins_code")
        if "ins_code" in binder.get_annotation_categories()
        else np.full(binder.array_length(), "")
    )
    cdr_mask = np.zeros(binder.array_length(), dtype=bool)
    print("--- CDR selection (source antibody numbering) ---")
    lengths = []
    for name, (first, last) in CDR_RANGES.items():
        m = (binder.res_id >= first) & (binder.res_id <= last)
        cdr_mask |= m
        labels = sorted(
            {(int(r), str(c)) for r, c in zip(binder.res_id[m], ins[m], strict=True)}
        )
        lengths.append(len(labels))
        shown = ",".join(f"{r}{c.strip()}" for r, c in labels)
        print(f"  {name} {first}-{last}: {len(labels)} residues  [{shown}]")
    print(f"  lengths {tuple(lengths)}; the 1G8 config cites {EXPECTED_CDR_LENGTHS}")
    if tuple(lengths) != EXPECTED_CDR_LENGTHS:
        print(
            "  NOTE lengths differ from the cited triple. Only refolded_apo_holo_cdr_rmsd\n"
            "       depends on this; refolded_apo_holo_rmsd does not. Adjust CDR_RANGES if\n"
            "       you want the CDR metric to match the campaign's own boundaries."
        )
    if cdr_mask.sum() == 0:
        raise ValueError("CDR selection is empty; CDR_RANGES does not match this numbering")

    # --- disulfide, before any renumbering changes the lookup -----------------
    a, b = BINDER_DISULFIDE
    ia = np.flatnonzero((binder.res_id == a) & (binder.atom_name == "SG"))
    ib = np.flatnonzero((binder.res_id == b) & (binder.atom_name == "SG"))
    if len(ia) == 1 and len(ib) == 1:
        d = float(np.linalg.norm(binder.coord[ia[0]] - binder.coord[ib[0]]))
        if d >= 2.5:
            raise ValueError(f"CYS{a}-CYS{b} SG-SG is {d:.2f} A; refusing to bond")
        binder.bonds.add_bond(int(ia[0]), int(ib[0]), BondType.SINGLE)
        print(f"  added disulfide CYS{a} SG - CYS{b} SG ({d:.2f} A)")
    else:
        print(f"  WARNING: CYS{a}/CYS{b} SG lookup failed; no disulfide declared")

    # --- renumber 1..N, clear insertion codes ---------------------------------
    # Insertion codes collapse onto their parent residue once dropped (52 and 52A both
    # become 52), which the CIF writer rejects as ambiguous. A residue-start scan treats an
    # insertion-code change as a boundary, so the two stay distinct and become consecutive.
    # Nothing downstream reads the binder's numbering: the oracle folds it from sequence and
    # every mask that touches it addresses the chain, not residue ids.
    starts = get_residue_starts(binder, add_exclusive_stop=True)
    new_ids = np.zeros(binder.array_length(), dtype=int)
    for i in range(len(starts) - 1):
        new_ids[starts[i] : starts[i + 1]] = i + 1
    binder.res_id = new_ids
    if "ins_code" in binder.get_annotation_categories():
        binder.set_annotation("ins_code", np.full(binder.array_length(), "", dtype="U1"))
    binder.chain_id = np.full(binder.array_length(), BINDER_CHAIN)

    # --- conditioning: framework sequence known, CDR sequence not -------------
    C_SEQ.set_annotation_from_ground_truth(binder, ~cdr_mask)
    C_IDX.set_annotation(binder, np.ones(binder.array_length(), dtype=bool))

    print(
        f"  binder: {binder.array_length()} atoms, {len(np.unique(binder.res_id))} residues, "
        f"C_SEQ on {int((~cdr_mask).sum())} atoms / absent on {int(cdr_mask.sum())}"
    )
    return binder


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading target via build_annotated_cifs.load_target() -- identical to the designs")
    target = load_target()
    print(f"    {target.array_length()} atoms, {len(np.unique(target.res_id))} residues\n")

    binder = load_native_binder()

    combined = combine(target.copy(), binder)
    print(f"\ncombined: {combined.array_length()} atoms, chains {list(np.unique(combined.chain_id))}")

    # Reported, not enforced: this is an rfd4-evaluate input, and the generation-input
    # validators expect coordinate conditioning on the binder that a control does not need.
    try:
        validate(combined, "pre-save")
    except ValueError as e:
        print(f"    note: generation-input validation reported:\n      {e}")

    to_cif_file(
        combined,
        str(OUT_CIF),
        save_standard_annotations=True,
        chain_disambiguation="chain_iid",
    )
    reloaded = as_atom_array_plus(parse(str(OUT_CIF), config=ANNOTATED_CIF)["assemblies"]["1"][0])
    compare_roundtrip(combined, reloaded, "native control")

    print("\n--- conditions on the saved file ---")
    for line in summarize_conditions(reloaded).splitlines():
        if line.strip():
            print("  " + line)

    print(f"\nwrote {OUT_CIF.relative_to(ROOT)}")
    print(
        "\nScore it with the campaign's own metric definitions:\n"
        "  PROJECT_ROOT=$PWD rfd4-evaluate --config-path $PWD/CSS/configs \\\n"
        "      --config-name css_native_control_evaluate \\\n"
        "      out_dir=./runs/native_control_annotated\n"
        "\nRead from ./runs/native_control_annotated/evaluation_metrics/per_example_metrics.csv:\n"
        "  oracle_ppi.min_ipae_coi            vs the designs' 15.8-24.5\n"
        "  refolded_apo_holo_rmsd.apo.rmsd    vs the 2.0 A criterion (designs: 4/12 passed)\n"
        "  refolded_apo_holo_cdr_rmsd.apo.rmsd\n"
        "  apo_holo_consistency.pass / .all_metrics_present\n"
        "\nall_metrics_present must be 1. A 0 means the control is indeterminate and the\n"
        "comparison is void -- the same trap the campaign's own smoke audit exists to catch."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
