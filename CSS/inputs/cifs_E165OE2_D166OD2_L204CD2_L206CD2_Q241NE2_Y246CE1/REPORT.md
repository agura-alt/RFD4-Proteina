# Annotated CIF build — "high constraint" arm

Built 2026-11-21. Hotspots **A/165 GLU OE2, A/166 ASP OD2, A/204 LEU CD2, A/206 LEU CD2,
A/241 GLN NE2, A/246 TYR CE1** (six).

Template: `CSS/inputs/cifs_L204CD2_L206CD2_Y246CE1/css_VHH_vobarilizumab_cdr123_L204CD2_L206CD2_Y246CE1.cif`
(each framework is diffed against **its own** template in that directory, not against vobarilizumab's).
Built by `CSS/build_annotated_cifs_highconstraint.py`, which reuses `CSS/build_annotated_cifs.py`
verbatim — the same pattern `build_annotated_cifs_groove.py` uses.

| File | Framework | Atoms | Bonds |
| --- | --- | --- | --- |
| `css_VHH_caplacizumab_cdr123_E165OE2_D166OD2_L204CD2_L206CD2_Q241NE2_Y246CE1.cif` | VHH_caplacizumab | 1781 | 1821 |
| `css_VHH_ozoralizumab_cdr123_E165OE2_D166OD2_L204CD2_L206CD2_Q241NE2_Y246CE1.cif` | VHH_ozoralizumab | 1776 | 1816 |
| `css_VHH_vobarilizumab_cdr123_E165OE2_D166OD2_L204CD2_L206CD2_Q241NE2_Y246CE1.cif` | VHH_vobarilizumab | 1778 | 1818 |
| `css_h-NbBCII10_cdr123_E165OE2_D166OD2_L204CD2_L206CD2_Q241NE2_Y246CE1.cif` | h-NbBCII10 | 1770 | 1809 |

## What changed relative to the template

Exactly one condition column: `condition_hotspot_1_atom` (3 atoms → 6 atoms). Verified
programmatically on all four files — atom identity, every coordinate, every bond, every
formal charge, and the other nine condition/annotation columns are **identical** to the
template:

```
condition columns identical: annotation_expsegmax_1_residue, annotation_expsegmin_1_residue,
  condition_CoordinateGroup_1_atom, condition_InternalCoordinate_1_atom,
  condition_cath-label_1_residue, condition_coordinate_1_atom,
  condition_index-group_1_token, condition_index_1_residue, condition_sequence_1_residue
condition columns changed  : condition_hotspot_1_atom
```

So `C_CRDGRP` (target group 0 = 560 **absolute** backbone atoms; binder group 1 = 332
internal-coordinate atoms), `C_SEQ` (all 140 target residues known; framework known, 3 CDR
sentinels free), `C_CATH`, the CDR windows (13-13, 9-10, 12-21) and both added disulfides
(CYS170–195, CYS226–264) are carried through unchanged. Chain H dock is still unconditioned
(`C_CRD` all-NaN on chain H).

## Atom selection

Per-atom Shrake–Rupley SASA on the isolated target (`biotite`, `vdw_radii="Single"`, 1000
points), sidechain atoms only:

| Residue | Chosen | SASA Å² | Alternatives (Å²) | Provenance of "same atom as the previous run" |
| --- | --- | --- | --- | --- |
| GLU 165 | **OE2** | 51.43 | OE1 41.47 · CG 14.37 · CD 7.61 · CB 7.00 | groove arm `cifs_E165OE2_L204CD2_L206CD2_Q241NE2`, user-confirmed 2026-11-19 |
| ASP 166 | **OD2** | 32.89 | CB 19.44 · OD1 12.00 · CG 5.31 | 1G8 positive control (`…_D166OD2_…`) |
| LEU 204 | **CD2** | 27.29 | CB 7.97 · CD1 0.12 · CG 0.00 | the named template run |
| LEU 206 | **CD2** | 21.25 | CD1 3.26 · all others 0.00 | the named template run |
| GLN 241 | **NE2** | 24.72 | **OE1 30.43** · CG 12.32 · CB 10.02 · CD 4.59 | 1G8 control + groove arm |
| TYR 246 | **CE1** | 17.87 | **OH 44.25** · CD1 3.74 · CE2 0.97 · CZ 0.60 | the named template run (the 1G8 control used OH; **CE1 confirmed 2026-11-21** — flag 2) |

All six are genuine epitope residues in `CSS_epitope_residue_report.csv`: 165 (98.5 % buried),
204 (97.7 %), 206 (94.9 %), 241 (92.3 %), 246 (91.8 %), 166 (82.8 %). That is epitope ranks
1, 2, 3, 4, 5 and 8 out of 25 — this arm takes essentially the whole top of the list.

## Validation performed

- `run_pipeline_validation`, `check_index_group_res_id_gaps`,
  `check_index_group_molecule_consistency` — passed pre-save and after reload, all four files.
- Round trip: 10 condition columns, all bonds, all charges identical, all four files.
- `CSS/verify_annotated_cifs.py` (re-pointed at this directory and hotspot set, run as a
  throwaway copy; the committed script still targets the original arm): **all checks passed
  for all four files** — chain A atom-identical to `target_crop.pdb` with both disulfides,
  chain H atom-identical to the repo example with its conditioning verbatim, coordinate
  groups 0/560 and 1/332 exactly as intended, hotspots exactly the six requested target atoms.
- `check_annotated_cif.sh`: 1 error, 3 warnings. Running it on the **template** as a control
  gives the *same* error and the *same three* warnings, including the same 8 `_struct_conn`
  records. Nothing new was introduced. The error is the known Windows-only
  `AddChiralityFeature` int32 bug (`chirality_feature.py:171`); the build script's scoped
  shim is what lets `run_pipeline_validation` actually run.

Environment note: the repo root has no `.venv`. These builds ran on
`CSS/.venv-fileprep/Scripts/python.exe`, the only environment on this machine carrying
`rfproteina` + `atomworks`. `check_annotated_cif.sh` needs `RFPROTEINA_PYTHON` pointed at it.

---

# Flags

## 1. Six hotspots at `min: 1.0` is a much harder gate than either predecessor

This is the main thing to decide before spending designs. `hotspot_obedience` is an
all-heavy-atom 4.5 Å nearest-neighbour criterion, and the campaign configs gate it at
`min: 1.0` — **every** hotspot contacted. That is now 6-of-6 rather than 3-of-3.

Geometry measured on the built files:

| | template arm (204/206/246) | groove arm (165/204/206/241) | **this arm (all six)** |
| --- | --- | --- | --- |
| hotspots | 3 | 4 | **6** |
| max pairwise span | 13.77 Å | 8.49 Å | **19.35 Å** |
| min pairwise span | 6.46 Å | 6.46 Å | 6.46 Å |
| max angular spread from centroid | 57.7° | 28.4° | **66.8°** |

The widest pair is **ASP166 OD2 ↔ TYR246 CE1 at 19.35 Å**, with GLN241 NE2 ↔ ASP166 OD2 at
15.05 Å close behind. A VHH paratope reaches roughly 20–25 Å, so the set is *co-reachable in
principle* and nothing here is structurally impossible — but it sits at the edge of that range
rather than comfortably inside it, and all six must land simultaneously. The template arm's own
config comment already warned to "expect a lower hotspot pass rate" at 13.77 Å and 57.7°; this
arm is wider on both numbers.

Burial proxy (heavy atoms within 8 Å): 165 OE2 **26**, 166 OD2 **44**, 204 CD2 77, 206 CD2 83,
241 NE2 48, 246 CE1 66. The two new additions are the two most solvent-exposed, least buried
atoms in the set — they are easy to touch individually, but they are what pulls the span to
19.35 Å.

**Recommendation, not acted on:** either relax `hotspot_adherence` to `min: 0.834` (5 of 6) or
`0.667` (4 of 6) in this arm's generate config, or accept that the headline pass rate will be
dominated by the gate rather than by design quality. Whichever you pick, score per-hotspot
offline — a pooled obedience average hides *which* atom a failing design dropped, and with six
hotspots that is the only actionable number. Do not compare this arm's `hotspot_adherence`
pass rate to the 3-hotspot arm's without saying so.

## 2. TYR 246: two different atoms in two previous runs — RESOLVED, CE1 confirmed

- The run named as the template (`cifs_L204CD2_L206CD2_Y246CE1`) used **CE1**.
- The 1G8 positive control (`…_D166OD2_Q241NE2_Y246OH`) used **OH**.

Both are "a previous run". **CE1** was used, because the user named that run as the template and
it is also the run supplying 204 and 206.

**Confirmed 2026-11-21: CE1 is intended. No change made; the files as built stand.**

Note for interpretation: CE1 is the *less* exposed of the two by a wide margin —
**OH 44.25 Å² vs CE1 17.87 Å²** — so this arm's TYR 246 contact is not the 1G8 control's, and
results should not be read as holding that aromatic hotspot constant against it. CE1 also sits
on the far side of the ring from OH, which is what makes the D166 OD2 ↔ Y246 CE1 pair the
19.35 Å widest span in flag 1. Measured both ways on `target_crop.pdb`, OH would have been the
tighter set — max span **17.36 Å** and max angular spread **58.3°**, against 19.35 Å and 66.8°
for CE1 (min span 6.46 Å either way). So the CE1 choice makes flag 1's gate about 2 Å harder.

Same shape of conflict, already settled in this arm's favour: **GLN 241** and **ASP 166** come
*only* from the 1G8 control, not from the named template. So this build already mixes two runs'
atom choices — flagged in the groove arm's report as its flag 2, and unchanged here.

## 3. "High constraint" collides with an existing "high" in this campaign

`CSS/configs/css_cdr123_L204CD2_L206CD2_Y246CE1_high_generate.yaml` and
`runs/css_cdr123_L204CD2_L206CD2_Y246CE1_high/` already exist, and there "high" means the
**high recycling arm** (300 trunk passes vs 180), not a hotspot count. Your "high constraint"
is a different axis entirely. I did **not** use "high" in any path — the directory and file
names encode the hotspot set the way every other arm here does
(`cifs_E165OE2_D166OD2_L204CD2_L206CD2_Q241NE2_Y246CE1`). Worth keeping that separation in the
generate config and run name too, or the two "high" arms will be indistinguishable in
`RUNS.md` and in any aggregation keyed on run name.

## 4. GLN 241 OE1 is more exposed than NE2 (30.43 vs 24.72 Å²)

Carried from the groove arm unchanged. NE2 was used because it is the 1G8 run's atom, but the
NE2/OE1 assignment in a glutamine amide is routinely flipped by refinement — crystallography
usually cannot distinguish N from O at this resolution. Both clear the 4.5 Å all-heavy-atom
criterion comfortably, so this is low-risk either way; if you want the single most exposed
atom, it is OE1.

## 5. E165 and D166 are now both in the set, which removes an earlier ambiguity

The groove arm's flag 3 recorded a worry that "165" might have been a slip for the 1G8
control's D166. That is moot here: this request names **both**, explicitly, as separate
residues. No inference was needed. Note for interpretation that they are adjacent residues
whose chosen atoms sit 10.28 Å apart, so they are not redundant, but they do pull in the same
direction.

## 6. Carried forward unchanged from the two previous builds, and still true

- **The lab's original 1G8 template CIF was never read** — it lives on a cluster path not
  mounted here. Conditioning is reconstructed from the repo examples.
- **`C_CRDGRP` "backbone-only" keeps the template's reading**: the *target* carries **absolute**
  backbone coordinates, pose pinned. `AGENTS.md` and
  `docs/how-to/input_prep/examples/condition_spec_ppi.yaml` give the PPI default as target
  backbone-only *internal* coordinates, pose free. You asked for the same `C_CRDGRP` as the
  template, so the template's reading is what was built — but the two differ in whether the
  antigen's global pose is pinned.
- **The two disulfides (CYS170–195, CYS226–264) are added, not preserved.** `target_crop.pdb`
  declares no SSBOND/CONECT and AtomWorks does not infer them. A deliberate chemistry change
  relative to the raw input, measured at 2.04 Å before bonding.
- **CDR lengths sample, they are not pinned** (13-13 / 9-10 / 12-21), so this is not comparable
  to the length-pinned 1G8 control.
- **`examples/antibody_design/vhh.yaml` would discard these hotspots** via
  `SampleNearbyHotspots`; pass `generation.dataset.sampling_transforms=null` or use a
  campaign-style config.
- **Chain H sits in the PMVK pose** and is not docked onto the CSS target in a viewer. `C_CRD`
  is NaN across chain H, so placement is sampled from noise.

## 7. Generate config — written 2026-11-21

`CSS/configs/css_6hot_E165OE2_D166OD2_L204CD2_L206CD2_Q241NE2_Y246CE1_generate.yaml`,
referenced on **`css_groove_E165OE2_L204CD2_L206CD2_Q241NE2_generate.yaml`** — the direct
predecessor: same four frameworks supplied as a directory, same medium/front_loaded recycling,
same metrics and pass criteria, and it already carries the 1000-design scale settings
(`compress_outputs: true`, `gen_njobs: 10`) that the arm before it did not.

Composed with Hydra and diffed leaf-by-leaf against that reference: **4 of 333 composed keys
differ** — `inputs`, `generation.val_set_name_override` (derived from `inputs`), `out_dir`,
`run_name` — and nothing else. Metrics, tags, pass criteria, recycling preset, `batch_size`,
`num_replicates`, `compress_outputs` and `gen_njobs` are untouched.

`CSS/preflight_metric_selection.py` on the new config returns output **identical** to the
reference: 60 steps / 180 trunk passes as expected, the same 11 inline metrics, 5 refolding
metrics, 2 oracles and 5 pass criteria on, and all five campaign-critical metrics `ok`.

**Design count: `num_replicates: 250` = 250 per framework = 1000 total**, matching the groove
arm so the 6-vs-4 hotspot comparison is at equal depth. *This was adopted from the predecessor,
not requested for this arm* — change `num_replicates` if a different depth was intended.

**`hotspot_adherence` left at `min: 1.0`**, i.e. 6-of-6, deliberately. Relaxing it buys nothing
that cannot be had offline: `hotspot_obedience.hotspot_obedience` is recorded as a raw fraction
on every design regardless of the gate, so a 5-of-6 or 4-of-6 analysis is a filter over the
output table after the fact, while changing the threshold here would make this the only arm in
the campaign whose hotspot gate differs. Expect a low pass rate on this criterion (flag 1).

**Name:** `css_6hot_` prefix, not `high` — see flag 3.

**Before submitting:** `gen_njobs: 10` is honoured by `scripts/launch_generate.sh`, which sizes
the SLURM array from the resolved config, but **not** by `CSS/submit_generate.sh`, which takes
`NJOBS` from the environment (default 1) and passes it as an explicit `gen_njobs=` override
that beats this file. With that launcher, set `NJOBS=10` and `--array=0-9` together — and note
that `CSS/submit_generate.sh` now refuses an array/`NJOBS` mismatch outright rather than
silently generating the wrong shard set (commit `153cf8c`). On Anvil,
`CSS/submit_generate_anvil.sh` already defaults `NJOBS` to 10 for a production submission.

(The groove report's §7 note that `CSS/submit_generate.sh` carries unreviewed uncommitted
modifications is now stale — those changes were committed in `153cf8c`.)
