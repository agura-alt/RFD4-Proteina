# Annotated CIF build — "wall only hydrophilic" arm

Built 2026-11-21. Hotspots **A/165 GLU OE2, A/166 ASP OD2, A/239 HIS CE1, A/241 GLN NE2**.

Template: `CSS/inputs/cifs_L204CD2_L206CD2_Y246CE1/css_VHH_vobarilizumab_cdr123_L204CD2_L206CD2_Y246CE1.cif`
Built by `CSS/build_annotated_cifs_wall.py` (reuses `CSS/build_annotated_cifs.py` and
`CSS/build_annotated_cifs_groove.py` verbatim).

> **HIS 239 atom switched NE2 → CE1 on instruction, 2026-11-21.** The first build of this
> arm used NE2; that directory (`cifs_E165OE2_D166OD2_H239NE2_Q241NE2/`) was deleted on
> instruction the same day, so this is the only "wall only hydrophilic" arm on disk.
> It can be regenerated from git history or by reverting one line in the build script.
> See flag 2.

| File | Framework | Atoms |
| --- | --- | --- |
| `css_VHH_caplacizumab_cdr123_E165OE2_D166OD2_H239CE1_Q241NE2.cif` | VHH_caplacizumab | 1781 |
| `css_VHH_ozoralizumab_cdr123_E165OE2_D166OD2_H239CE1_Q241NE2.cif` | VHH_ozoralizumab | 1776 |
| `css_VHH_vobarilizumab_cdr123_E165OE2_D166OD2_H239CE1_Q241NE2.cif` | VHH_vobarilizumab | 1778 |
| `css_h-NbBCII10_cdr123_E165OE2_D166OD2_H239CE1_Q241NE2.cif` | h-NbBCII10 | 1770 |

As in the previous arm, all four are diffed against their **own** framework's template
in `cifs_L204CD2_L206CD2_Y246CE1/`; in every case only `condition_hotspot_1_atom` differs.

## What changed relative to the template

Exactly one condition column: `condition_hotspot_1_atom` (3 atoms → 4 atoms).
Verified programmatically per framework — atom identity, every coordinate, every bond
(1821 / 1816 / 1818 / 1809), all formal charges, and the other nine condition/annotation
columns are **identical** to that framework's template:

```
condition columns identical: annotation_expsegmax_1_residue, annotation_expsegmin_1_residue,
  condition_CoordinateGroup_1_atom, condition_InternalCoordinate_1_atom,
  condition_cath-label_1_residue, condition_coordinate_1_atom,
  condition_index-group_1_token, condition_index_1_residue, condition_sequence_1_residue
condition columns changed  : condition_hotspot_1_atom
```

So `C_CRDGRP` (target group 0 = 560 absolute backbone atoms; binder group 1 = 332
internal-coordinate atoms), `C_SEQ` (all 140 target residues known; framework known,
3 CDR sentinels free), `C_CATH`, the CDR windows (13-13, 9-10, 12-21) and both added
disulfides are carried through unchanged. Bonds and charges are preserved by
construction, not re-derived.

## Atom selection

Per-atom Shrake–Rupley SASA on the isolated target (`biotite`, `vdw_radii="Single"`,
1000 points). Chosen atom in bold.

| Residue | Chosen | SASA Å² | Alternatives (Å²) | Provenance |
| --- | --- | --- | --- | --- |
| GLU 165 | **OE2** | 51.43 | OE1 41.47 · CG 14.37 · CD 7.61 · CB 7.00 | previous run `cifs_E165OE2_...` |
| ASP 166 | **OD2** | 32.89 | CB 19.44 · OD1 12.00 · CG 5.31 | 1G8 positive control (`D166OD2`) |
| HIS 239 | **CE1** | 26.93 | CD2 24.39 · NE2 18.26 · **ND1 2.52** | **no precedent — user's choice, flag 2** |
| GLN 241 | **NE2** | 24.72 | OE1 30.43 · CG 12.32 · CB 10.02 · CD 4.59 | 1G8 control + `cifs_E165OE2_...` |

Three of the four atoms are literal carry-overs, as you asked. Only HIS 239 required a
choice, because no previous run has ever used residue 239. CE1 is the most
solvent-exposed atom of the residue.

All four residues are genuine epitope residues in `CSS_epitope_residue_report.csv`, and
all four are flagged `hydrophobic = False` there. Burial ranks among the 25 epitope
residues: 165 (98.5 %, rank 1), 239 (95.7 %, rank 3), 241 (92.3 %, rank 4), 166
(82.8 %, rank 9). The two leucines of the previous two arms — the only hydrophobic
members — are both dropped. Note that the *residues* are all hydrophilic but the chosen
HIS 239 **atom** is now a ring carbon; see flag 1.

## Validation performed

- `run_pipeline_validation`, `check_index_group_res_id_gaps`,
  `check_index_group_molecule_consistency` — passed pre-save **and** after reload, on
  all four files.
- Round trip: 10 condition columns, all bonds, all charges identical.
- Independent read-back through `ANNOTATED_CIF` confirms 4 hotspot atoms per file,
  exactly `A/165/GLU/OE2, A/166/ASP/OD2, A/239/HIS/CE1, A/241/GLN/NE2`.
- `check_annotated_cif.sh`: **1 error, 3 warnings** on each of the four. Running the
  checker on the **template** as a control gives the same error and the same three
  warnings, including the same 8 `_struct_conn` records and the same
  `+1: 24 / -1: 17` charge inventory. Nothing new was introduced. The error is the
  known Windows-only `AddChiralityFeature` int32 bug
  (`src/rfproteina/transforms/chirality_feature.py:171`); the build script's scoped shim
  is what lets `run_pipeline_validation` actually run.

---

# Flags

## 1. With CE1, the arm name no longer describes one of its four contacts

All four hotspot *residues* are non-hydrophobic, which is what made the name "wall only
hydrophilic" accurate. But the atom now nominated on HIS 239 is **CE1, a ring carbon**,
not a nitrogen. `hotspot_obedience` is an all-heavy-atom 4.5 Å nearest-neighbour
criterion, so it does not care either way — but if the arm's purpose is to test an
all-polar epitope against the leucine-containing arms, this contact is no longer
carrying that hypothesis. Three of four contacts are polar (OE2, OD2, NE2); the fourth
is aromatic carbon.

This is a deliberate choice you made (raw exposure over H-bond capability) and it is
defensible — CE1 is the most exposed atom on the residue at 26.93 Å², and histidine's
ring face is a real binding determinant. Recording it because the arm name and the arm
contents now disagree slightly, and six months from now that will not be obvious.

**Geometry is essentially unchanged by the switch.** All numbers recomputed with one
consistent method (angular spread about the target centroid; burial proxy = heavy atoms
within 8 Å):

| | arm 1 `204/206/246` | arm 2 `165/204/206/241` | **this arm (CE1)** | first draft (NE2) | 1G8 control |
| --- | --- | --- | --- | --- | --- |
| hotspots | 3 | 4 | **4** | 4 | 3 |
| max pairwise span | 13.77 Å | 8.49 Å | **15.05 Å** | 15.05 Å | 17.36 Å |
| min pairwise span | 6.46 Å | 6.46 Å | **4.80 Å** | 5.01 Å | 7.96 Å |
| angular spread | 57.7° | 28.4° | **48.1°** | 48.1° | 58.3° |
| burial proxy | 77 / 83 / 66 | 26 / 77 / 83 / 48 | **26 / 44 / 43 / 48** | 26 / 44 / 37 / 48 | 44 / 48 / 61 |

Every atom sits in the 26–48 neighbour range; arm 1 had two atoms at 77 and 83. So this
epitope is substantially **less recessed** than arm 1 and **wider** than arm 2, landing
close to the 1G8 positive control on span and angular spread. On pure reachability this
should be the easiest of the three CSS arms to satisfy.

**But chemistry and geometry pull in opposite directions here.** Removing LEU 204 and
LEU 206 removes the only buried apolar anchors in the set. Protein interfaces are driven
largely by buried apolar surface, and a near-all-polar hotspot set asks the model to pay
desolvation costs on these contacts without that compensation. 165 and 166 are both
acidic, so designs will be pushed toward an Arg/Lys-rich paratope — a common failure mode
for *specificity*, which often refolds with lower confidence even when
`hotspot_obedience` passes. **Expect `hotspot_obedience` to pass more readily here than
in arm 1, while iPAE / pAE_interaction may not improve, and may get worse.** Do not read
a high hotspot pass rate on this arm as a better binder. If the intent was a hydrophilic
rim around a retained hydrophobic core, keeping one leucine (204 CD2, 27.3 Å²) would
test that directly.

## 2. HIS 239 remains the one atom with no precedent

No previous run used residue 239, so "the same atom as the previous run" has no referent
here. CE1 was selected by you on 2026-11-21 as the most exposed atom (26.93 Å²), in
preference to NE2 (18.26 Å², the most exposed nitrogen), which this arm's first build
used.

- **ND1 was never a candidate: it is buried at 2.52 Å²**, the exact HIS analogue of the
  LEU CD1 problem resolved in the first build. Under the 4-of-4 `hotspot_obedience` gate
  it would have been close to structurally unreachable.
- **CD2 (24.39 Å²) is the near-tie with CE1** and sits on the opposite side of the
  imidazole. Either would serve if you ever want the other ring face.
- **The superseded NE2 build was deleted** on instruction, 2026-11-21, having been
  complete and separately validated. The reason was collision risk, not a defect:
  campaign configs take `inputs` as a *directory*, and two directories differing by a
  single atom are easy to mis-point at submission time. To recreate it, set HIS 239 back
  to `NE2` in `HOTSPOTS` and `TAG`/`OUT_DIR` in `CSS/build_annotated_cifs_wall.py`.
- HIS protonation/tautomer state is not modelled here. This no longer matters for the
  hotspot itself (CE1 is a carbon, unaffected by the ND1/NE2 assignment ambiguity) —
  which is a minor side benefit of the switch.

## 3. The requested hotspot list is out of numeric order — confirm no transcription slip

You wrote them as **165 E, 166 D, 241 Q, 239 H** — 241 before 239. `C_HOT` is an
unordered atom mask, so ordering has **no effect whatsoever** on the build, and the
files are written in numeric order. I raise it only because an out-of-order last pair is
the shape a copy/paste slip takes, and 239/241 are 4.80 Å apart — the tightest pair in
the set. Residue identities all check out against the structure (165 GLU, 166 ASP,
239 HIS, 241 GLN), so if the letters are right the numbers almost certainly are too.
No action taken.

## 4. GLN 241 OE1 is still more exposed than NE2 (30.43 vs 24.72 Å²)

Carried forward unchanged from the previous arm's flag 4. NE2 is used because it is the
1G8 control's atom, which is what you asked for. The NE2/OE1 assignment in a glutamine
amide is routinely flipped by refinement. Low risk either way; OE1 is the more exposed
atom if you ever want to switch. Note that having now switched HIS 239 to its
most-exposed atom, GLN 241 is the remaining hotspot *not* on its most-exposed atom — if
the principle behind the CE1 switch was "take maximum exposure", consistency would argue
for OE1 here too.

## 5. Four hotspots again — `hotspot_obedience >= 1.0` is a 4-of-4 gate

The template config gates at `min: 1.0`, i.e. *all* hotspots contacted. As with the
previous arm this is a strictly harder gate than arm 1's 3-of-3. Here it runs in the
*same* direction as the geometry rather than against it (flag 1: the epitope is more
exposed), but it still means the three CSS arms are not gated comparably. If you want
like-for-like pass rates across arms, either hold the hotspot count fixed or relax
`min` to 0.75.

## 6. E165 and D166 are sequence-adjacent but not geometrically redundant

Worth recording, since spending two of four hotspots on consecutive residues looks like
it wastes half the budget. It does not: their chosen sidechain atoms point in different
directions and are **10.28 Å** apart — wider than the 239–241 pair (4.80 Å). All six
pairwise distances: 165–166 10.28, 165–239 9.26, 165–241 8.49, 166–239 12.18,
166–241 15.05, 239–241 4.80 Å. The 15.05 Å maximum is comfortably inside a VHH
paratope's ~20–25 Å reach.

## 7. Carried forward from the previous builds, unchanged and still true

- **The lab's original 1G8 template CIF was never read** — it lives on a cluster path
  not mounted here. Conditioning is reconstructed from the repo examples, and this arm
  inherits it transitively through `cifs_L204CD2_L206CD2_Y246CE1`. (Arm 1 flag 5.)
- **`C_CRDGRP` "backbone-only" keeps the template's reading**: the *target* carries
  **absolute** backbone coordinates, pose pinned. `AGENTS.md` and
  `docs/how-to/input_prep/examples/condition_spec_ppi.yaml` state the PPI default as
  target backbone-only *internal* coordinates, pose free. You asked for the same
  `C_CRDGRP` as the template, so the template's reading is what was built — but the two
  differ in whether the antigen's global pose is pinned. (Arm 1 flag 6, arm 2 flag 6.)
- **The two disulfides (CYS170–195, CYS226–264) are added, not preserved** —
  `target_crop.pdb` has no SSBOND/CONECT and AtomWorks does not infer them from
  distance. "Preserve bonds" therefore means preserving the *template's* bond set, which
  is what was verified; relative to the raw PDB this remains a deliberate chemistry
  addition. (Arm 1 flag 7.)
- **CDR lengths sample, they are not pinned** (13-13 / 9-10 / 12-21), so this is not
  comparable to the length-pinned 1G8 control. (Arm 1 flag 9.)
- **`examples/antibody_design/vhh.yaml` would discard these hotspots** via
  `SampleNearbyHotspots`; pass `generation.dataset.sampling_transforms=null` or use a
  campaign-style config. (Arm 1 flag 11.)
- **Chain H sits in the PMVK pose** and is not docked onto the CSS target in a viewer.
  `C_CRD` is NaN across chain H, so placement is sampled from noise. (Arm 1 flag 13.)
- **The Windows `chirality_feature.py:171` bug** (`np.random.randint(0, 2**32)` under
  the legacy int32 default) still blocks `run_pipeline_validation` unshimmed, and still
  deserves a one-line PR (`np.random.default_rng()`). (Arm 1 flag 10.)

## 8. Generate config — built 2026-11-21

`CSS/configs/css_wallphilic_E165OE2_D166OD2_H239CE1_Q241NE2_generate.yaml`, referenced on
**`css_6hot_E165OE2_D166OD2_L204CD2_L206CD2_Q241NE2_Y246CE1_generate.yaml`** — the most
recent predecessor and the one that had already settled the 1000-design scale settings
(`compress_outputs: true`, `gen_njobs: 10`), so referencing it keeps the delta list to
three keys instead of five. The groove arm would have served equally well, but its header
carries a stale contradiction (its "Design count" section says `compress_outputs` and
`gen_njobs` are not set, while its body sets both); referencing the 6hot file avoids
propagating that text. **Worth fixing in the groove config independently.**

Verified by Hydra composition diff: **418 composed keys on each side, 4 differing** —
`inputs`, `out_dir`, `run_name`, and `generation.val_set_name_override`, which is derived
from `inputs` and tracks it automatically. Metrics, pass criteria, recycling preset,
`num_replicates`, `batch_size` and both scale settings are untouched.

`CSS/preflight_metric_selection.py` passes on the new config: 60 sampling steps /
180 trunk passes as expected for PPI medium front-loaded, 11 inline metrics on,
5 refolding metrics on, both oracles on, 5 pass criteria on, and all five
campaign-critical metrics (`hotspot_obedience`, `interchain_interactions`,
`antibody_contacts`, `design_motif_allatom_rmsd`, `sequence_adherence`) resolve to ON.

**Design count: `num_replicates: 250` = 250 per framework = 1000 total**, since
`num_replicates` is per input and `inputs` is a directory of 4 CIFs. This matches the
groove and 6hot arms, so all three multi-framework arms sit at equal depth. *Assumption:*
250 was adopted from the predecessors, not requested for this arm.

**Before submitting:** the previous arm's flag 7 still applies. `CSS/submit_generate.sh`
**ignores `gen_njobs` from the config** and passes its own `NJOBS` (default 1) as an
explicit override that beats this file, so set `NJOBS=10` and `--array=0-9` together if
you use that launcher. `scripts/launch_generate.sh` reads it from the resolved config and
needs no such care. `CSS/submit_generate.sh` also still has uncommitted local
modifications that were not reviewed as part of this build.
