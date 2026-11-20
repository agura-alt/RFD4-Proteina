# Annotated CIF build — "spanning the groove" arm

Built 2026-11-19. Hotspots **A/165 GLU OE2, A/204 LEU CD2, A/206 LEU CD2, A/241 GLN NE2**.

Template: `CSS/inputs/cifs_L204CD2_L206CD2_Y246CE1/css_VHH_vobarilizumab_cdr123_L204CD2_L206CD2_Y246CE1.cif`
Built by `CSS/build_annotated_cifs_groove.py` (reuses `CSS/build_annotated_cifs.py` verbatim).

| File | Framework | Atoms |
| --- | --- | --- |
| `css_VHH_caplacizumab_cdr123_E165OE2_L204CD2_L206CD2_Q241NE2.cif` | VHH_caplacizumab | 1781 |
| `css_VHH_ozoralizumab_cdr123_E165OE2_L204CD2_L206CD2_Q241NE2.cif` | VHH_ozoralizumab | 1776 |
| `css_VHH_vobarilizumab_cdr123_E165OE2_L204CD2_L206CD2_Q241NE2.cif` | VHH_vobarilizumab | 1778 |
| `css_h-NbBCII10_cdr123_E165OE2_L204CD2_L206CD2_Q241NE2.cif` | h-NbBCII10 | 1770 |

All four are diffed against their **own** framework's template; in every case only
`condition_hotspot_1_atom` differs.

## What changed relative to the template

Exactly one condition column: `condition_hotspot_1_atom` (3 atoms → 4 atoms).
Verified programmatically — atom identity, all 1778 coordinates, all 1818 bonds, all
formal charges, and the other nine condition/annotation columns are **identical** to
the template:

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
disulfides are carried through unchanged.

## Atom selection

Per-atom Shrake–Rupley SASA on the isolated target (`biotite`, `vdw_radii="Single"`,
1000 points):

| Residue | Chosen | SASA Å² | Alternatives |
| --- | --- | --- | --- |
| GLU 165 | **OE2** | **51.43** | OE1 41.47 · CG 14.37 · CD 7.61 · CB 7.00 |
| LEU 204 | **CD2** | **27.29** | CB 7.97 · CD1 0.12 · CG 0.00 |
| LEU 206 | **CD2** | **21.25** | CD1 3.26 · all others 0.00 |
| GLN 241 | **NE2** | **24.72** | **OE1 30.43** · CG 12.32 · CB 10.02 · CD 4.59 |

Provenance of "the same atom used in the previous run":

- **204 CD2, 206 CD2** — literally the previous run (`cifs_L204CD2_L206CD2_Y246CE1`).
- **241 NE2** — from the 1G8 positive control
  `css_1g8pin_1G8_cdr123_L13-10-14_D166OD2_Q241NE2_Y246OH_generate.yaml`, a *different*
  previous run. See flag 2.
- **165 OE2** — no previous run used 165. OE2 is the Glu counterpart of the 1G8 run's
  D166 **OD2**, and is the most exposed atom of the residue. See flag 3.

All four residues are genuine epitope residues in `CSS_epitope_residue_report.csv` —
in fact the **top four by burial** out of 25: 165 (98.5 %), 204 (97.7 %), 206 (94.9 %),
241 (92.3 %). TYR 246 (91.8 %) was rank 5 and is the one dropped.

## Validation performed

- `run_pipeline_validation`, `check_index_group_res_id_gaps`,
  `check_index_group_molecule_consistency` — passed pre-save and after reload.
- Round trip: 10 condition columns, 1818 bonds, all charges identical.
- `check_annotated_cif.sh`: 1 error, 3 warnings. Running the checker on the **template**
  as a control gives the *same* error and the *same three* warnings, including the same
  8 `_struct_conn` records. Nothing new was introduced. The error is the known
  Windows-only `AddChiralityFeature` int32 bug (`chirality_feature.py:171`); the build
  script's scoped shim is what lets `run_pipeline_validation` actually run.

---

# Flags

## 1. The colloquial name points the wrong way — this set spans *less*, not more

"Spanning the groove" suggests a wider, more wrap-around epitope than the previous arm.
Measured on the final atom set, it is the opposite — the new tetrad is **tighter and
flatter** than the triad it replaces:

| | Previous arm (204 CD2 / 206 CD2 / 246 CE1) | **This arm** (165 OE2 / 204 CD2 / 206 CD2 / 241 NE2) |
| --- | --- | --- |
| max angular spread from centroid | 57.7° | **28.4°** |
| max pairwise span | 13.77 Å | **8.49 Å** |
| min pairwise span | 6.46 Å | 6.46 Å |
| burial proxy (heavy atoms within 8 Å) | 78 / 84 / 67 | 27 / **78 / 84** / 49 |

Adding GLU 165 (exposed, only 27 heavy neighbours) and GLN 241 (49) *flattens* the set:
both new atoms are polar, solvent-facing and sit on the same face, while the atom that
was pulling the old set around the antigen — TYR 246 CE1 — is the one being removed.
This should be an **easier** `hotspot_obedience` target than the previous arm, not a
harder one. If the intent behind the name was to force a binder across the groove,
this hotspot set does not do that, and you may want 246 kept in. Flagged, not acted on.

## 2. "The same atom as the previous run" resolves to two different runs

204 and 206 come from the `L204CD2_L206CD2_Y246CE1` run named as the template. 241 does
not appear in that run at all — `Q241 NE2` comes from the earlier 1G8 positive control.
I used NE2. Confirm that is the run you meant.

## 3. GLU 165 has no precedent in any previous run — RESOLVED, E165 OE2 confirmed

You asked for "the same atom that was used in the previous run", but no run has ever
used residue 165. OE2 is the defensible choice (most exposed; direct counterpart of the
1G8 run's D166 OD2), but it was **an inference, not a carry-over**.

The nearby alternative was the 1G8 template's acidic hotspot **ASP 166 OD2** — the
residue immediately after the 165 requested, also acidic, also in the epitope report
(82.8 % burial) — which would have made a 165/166 slip plausible.

**Confirmed 2026-11-19: E165 OE2 is intended, not a slip for D166. No change made; the
file as built stands.** Note for interpretation: E165 and D166 are different epitope
positions, so this arm's acidic contact is not the 1G8 control's, and results against
that control should not be read as holding the acidic hotspot constant.

## 4. GLN 241 OE1 is more exposed than NE2 (30.43 vs 24.72 Å²)

NE2 was used because it is the 1G8 run's atom. But the NE2/OE1 assignment in a glutamine
amide is frequently flipped by refinement — crystallography usually cannot distinguish
N from O at this resolution — exactly the ambiguity that forced the CD1→CD2 switch on
the leucines in the previous build. Both atoms clear the 4.5 Å all-heavy-atom
`hotspot_obedience` criterion comfortably, so this is low-risk either way, but if you
want the single most exposed atom, it is OE1.

## 5. Four hotspots changes what `hotspot_obedience >= 1.0` means

The template config gates at `min: 1.0`, i.e. *all* hotspots contacted. That is now a
**4-of-4** requirement rather than 3-of-3, which is a strictly harder gate even though
the epitope itself got easier (flag 1). The two effects run in opposite directions, so
do not read the pass rate against the previous arm without accounting for both. If you
want the gate comparable, either keep 3 hotspots or relax `min` to 0.75.

## 6. Carried forward from the previous build, unchanged and still true

- **The lab's original 1G8 template CIF was never read** — it lives on a cluster path
  not mounted here. Conditioning is reconstructed from the repo examples. (Prev. flag 5.)
- **`C_CRDGRP` "backbone-only" keeps the template's reading**: the *target* carries
  **absolute** backbone coordinates, pose pinned. `AGENTS.md` and
  `docs/how-to/input_prep/examples/condition_spec_ppi.yaml` state the PPI default as
  target backbone-only *internal* coordinates, pose free. You asked for "the same
  C_CRDGRP as the template", so the template's reading is what was built — but the two
  differ in whether the antigen's global pose is pinned. (Prev. flag 6.)
- **The two disulfides (CYS170–195, CYS226–264) are added, not preserved** —
  `target_crop.pdb` has no SSBOND/CONECT and AtomWorks does not infer them. This is a
  deliberate chemistry change relative to the raw input. (Prev. flag 7.)
- **CDR lengths sample, they are not pinned** (13-13 / 9-10 / 12-21), so this is not
  comparable to the length-pinned 1G8 control. (Prev. flag 9.)
- **`examples/antibody_design/vhh.yaml` would discard these hotspots** via
  `SampleNearbyHotspots`; pass `generation.dataset.sampling_transforms=null` or use a
  campaign-style config. (Prev. flag 11.)
- **Chain H sits in the PMVK pose** and is not docked onto the CSS target in a viewer.
  `C_CRD` is NaN across chain H, so placement is sampled from noise. (Prev. flag 13.)

## 7. Generate config — built 2026-11-19

`CSS/configs/css_groove_E165OE2_L204CD2_L206CD2_Q241NE2_generate.yaml`, referenced on
**`css_cdr123_L204CD2_L206CD2_Y246CE1_generate.yaml`** (the direct predecessor: same four
frameworks via directory input, same medium/front_loaded recycling, same metrics and
gates). Composed with Hydra and diffed leaf-by-leaf against that reference: **4 of 408
composed keys differ** — `inputs`, `out_dir`, `run_name`, `num_replicates` — and nothing
else. Metrics, pass criteria, recycling preset and `batch_size` are untouched.

**Design count: `num_replicates: 250` = 250 per framework = 1000 total.** `num_replicates`
is *per input* and `inputs` is a directory of 4 CIFs, so the total is 4 × 250. Confirmed
2026-11-19 (the first draft read "250 across the four frameworks" as 250 *total* and used
63; that was wrong and has been corrected).

### Scale settings, adopted from the other 1000-design run

On instruction (2026-11-19) this arm now matches `css_vobarilizumab_1000_generate.yaml`,
the campaign's other 1000-design run, on both scale settings:

| | this arm | 1000-design sibling | reference (`..._Y246CE1_generate`) |
| --- | --- | --- | --- |
| `compress_outputs` | **true** | true | false |
| `gen_njobs` | **10** | 10 | 1 |
| `num_replicates` | 250 (×4 inputs) | 1000 (×1 input) | 100 |
| `batch_size` | 2 | 2 | 2 |

So the composed config now differs from its reference in **6 of 408 keys**: the four
above minus `batch_size`, plus `inputs`, `out_dir` and `run_name`. Verified by Hydra
composition diff.

**`gen_njobs: 10` behaves differently under the two launchers — check before submitting.**

- `scripts/launch_generate.sh` reads `gen_njobs` from the resolved config when `--njobs`
  is not passed, and sizes the SLURM array from it. The config's 10 is honoured.
- `CSS/submit_generate.sh` does not read the config. It takes `NJOBS` from the
  environment, defaulting to **1**, and passes it to `rfd4-generate` as an explicit
  `gen_njobs=` override — which **beats** the value in this file. Submitting without
  setting it would silently run this arm at `gen_njobs=1` on a single array task.
  That script's own header documents the pattern for the 1000-design probes:
  set `NJOBS=10` and `--array=0-9` together, so the array size and `gen_njobs` agree.

Note also that `CSS/submit_generate.sh` has uncommitted local modifications that were
not reviewed as part of this build.
