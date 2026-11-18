"""
Rank target residues at a nanobody-target interface.

Step 1: enumerate contacts (heavy-atom distance cutoff, and 10 A CB cutoff).
Step 2: per-residue buried SASA = SASA(target alone) - SASA(target in complex).

Usage:
    python interface_hotspots.py complex.cif A B
        where A = target chain, B = nanobody chain
"""

import sys
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

CIF, TARGET_CHAIN, NB_CHAIN = sys.argv[1], sys.argv[2], sys.argv[3]
HEAVY_CUTOFF = 4.5   # A, any-heavy-atom contact
CB_CUTOFF = 10.0     # A, matches the RFdiffusion training definition

# --- load -------------------------------------------------------------
# use_author_fields=False (the default in biotite 2.x) gives label_seq_id /
# label_asym_id, which is the numbering RFD4 selection strings use.
# Set it True if you want to compare against what PyMOL shows you.
block = pdbx.CIFFile.read(CIF)
arr = pdbx.get_structure(block, model=1, use_author_fields=False)

arr = arr[struc.filter_amino_acids(arr)]          # drop waters, ions, ligands
arr = arr[arr.element != "H"]                     # heavy atoms only

tgt = arr.chain_id == TARGET_CHAIN
nb = arr.chain_id == NB_CHAIN
assert tgt.sum() and nb.sum(), "check your chain IDs"

# --- step 1: contacts -------------------------------------------------
cell = struc.CellList(arr[nb], cell_size=HEAVY_CUTOFF)
hits = cell.get_atoms(arr[tgt].coord, radius=HEAVY_CUTOFF)
atom_in_contact = (hits != -1).any(axis=1)

# collapse to residues: a residue is in contact if any heavy atom is
res_contact = struc.apply_residue_wise(arr[tgt], atom_in_contact, np.any)

# count contacts contributed by side-chain atoms only (beyond CB)
sidechain = ~np.isin(arr[tgt].atom_name, ["N", "CA", "C", "O", "CB"])
res_sc_contacts = struc.apply_residue_wise(
    arr[tgt], (atom_in_contact & sidechain).astype(int), np.sum
)
res_bb_contacts = struc.apply_residue_wise(
    arr[tgt], (atom_in_contact & ~sidechain).astype(int), np.sum
)

# the wider CB-based definition the model was trained on
cb = arr[tgt][np.isin(arr[tgt].atom_name, ["CB", "CA"])]   # CA stands in for GLY
cb = cb[~((cb.atom_name == "CA") & (cb.res_name != "GLY"))]
cell_nb_cb = struc.CellList(arr[nb], cell_size=CB_CUTOFF)
cb_hits = cell_nb_cb.get_atoms(cb.coord, radius=CB_CUTOFF)
cb_near = (cb_hits != -1).any(axis=1)
cb_lookup = dict(zip(cb.res_id, cb_near))

# --- step 2: buried SASA ---------------------------------------------
sasa_free = struc.sasa(arr[tgt], point_number=1000)           # target alone
sasa_bound_all = struc.sasa(arr, atom_filter=tgt, point_number=1000)
sasa_bound = sasa_bound_all[tgt]                              # NaN outside filter

free = struc.apply_residue_wise(arr[tgt], sasa_free, np.nansum)
bound = struc.apply_residue_wise(arr[tgt], sasa_bound, np.nansum)
buried = free - bound

res_ids, res_names = struc.get_residues(arr[tgt])

# --- report -----------------------------------------------------------
order = np.argsort(-buried)
print(f"{'res':>6} {'name':>5} {'buriedA2':>9} {'freeA2':>8} {'%buried':>8} "
      f"{'sc_cts':>7} {'bb_cts':>7} {'4.5A':>5} {'10A_CB':>7}")
for i in order:
    if buried[i] < 1.0 and not res_contact[i]:
        continue
    pct = 100 * buried[i] / free[i] if free[i] > 0 else 0.0
    print(f"{res_ids[i]:>6} {res_names[i]:>5} {buried[i]:>9.1f} {free[i]:>8.1f} "
          f"{pct:>8.1f} {res_sc_contacts[i]:>7} {res_bb_contacts[i]:>7} "
          f"{str(bool(res_contact[i])):>5} "
          f"{str(bool(cb_lookup.get(res_ids[i], False))):>7}")

print(f"\ntotal buried on target: {buried.sum():.0f} A^2")
print(f"contact residues (4.5 A heavy): {int(res_contact.sum())}")
print(f"residues within {CB_CUTOFF} A CB: {int(sum(cb_lookup.values()))}")
