# Annotated CIF build — CSS target × four repo VHH frameworks

Built 2026-11-16. Hotspots **A/204 LEU CD2, A/206 LEU CD2, A/246 TYR CE1**.

## Files

| File | Framework | Atoms |
| --- | --- | --- |
| `css_VHH_caplacizumab_cdr123_L204CD2_L206CD2_Y246CE1.cif` | VHH_caplacizumab | 1781 |
| `css_VHH_ozoralizumab_cdr123_L204CD2_L206CD2_Y246CE1.cif` | VHH_ozoralizumab | 1776 |
| `css_VHH_vobarilizumab_cdr123_L204CD2_L206CD2_Y246CE1.cif` | VHH_vobarilizumab | 1778 |
| `css_h-NbBCII10_cdr123_L204CD2_L206CD2_Y246CE1.cif` | h-NbBCII10 | 1770 |

Built by `CSS/build_annotated_cifs.py`, verified by `CSS/verify_annotated_cifs.py`.
Uncompressed `.cif` (the repo examples ship `.cif.zst`; both are accepted).

## Provenance

- **Target, chain A** — `CSS/target_crop.pdb`, parsed with the `DISTILLATION` preset
  (`load_standard_annotations=False`, since AtomWorks refuses that flag on PDB input).
  1134 atoms, 140 residues, auth numbering 151–290, no gaps, no insertion codes,
  no altlocs, no HETATM, no missing sidechain atoms. `add_missing_atoms=False`, so
  nothing was repaired or rebuilt.
- **Binder, chain H** — carried over verbatim from
  `examples/antibody_design/structures/vhh/pmvk_3CH4__<framework>_cdr123.cif.zst`.
  Chain H is atom-identical between the `pmvk_3CH4` and `mzb1_7AAH` copies of each
  framework; only the (unconditioned) dock pose differs, so `pmvk_3CH4` was taken
  arbitrarily.
- **Environment** — `CSS/.venv-fileprep` (Python 3.12.10, atomworks 3.0.4,
  biotite 1.6.0, rfd4-proteina 0.0.3, torch 2.15.0 CPU). File prep only; no GPU,
  no checkpoints, nothing runnable.
- No tracked repository file was modified. `git status` is unchanged (`?? CSS/`).

## Conditioning

Reproduces the scheme of the repo's own VHH example inputs and of the campaign
template `css_1g8pin_1G8_cdr123_L13-10-14_D166OD2_Q241NE2_Y246OH_generate.yaml`:

| Condition | Chain A (target) | Chain H (binder) |
| --- | --- | --- |
| `C_CRDGRP` | group **0** — 560 backbone atoms | group **1** — 332 backbone atoms |
| `C_CRD` (absolute) | finite on those 560 atoms | **NaN throughout** — dock not conditioned |
| `C_INTCRD` (internal) | none | finite on the 332 group-1 atoms |
| `C_SEQ` | all 140 residues known | framework known; 3 CDR sentinels free |
| `C_HOT` | 3 atoms (204 CD2, 206 CD2, 246 CE1) | none |
| `C_CATH` | `x.x.x.x` | `2.60.40.x` (Ig fold) |

CDR expandable segments carried over unchanged: **(13,13), (9,10), (12,21)**.
Realized binder length is therefore 83 + 34…44 = **117–127 residues**.

## Validation performed

For every file, both before saving and after reloading with `ANNOTATED_CIF`:

- `run_pipeline_validation`, `check_index_group_res_id_gaps`,
  `check_index_group_molecule_consistency` — all passed.
- Round trip compared atom-by-atom on atom identity: 10 condition columns, the full
  bond set, and all formal charges identical.
- Independent diff against both sources: target atom identity, coordinates and
  charges identical to `target_crop.pdb`; chain H identity, charges, and all six
  conditioning/segment columns identical to the repo example.
- `skills/rfd4-proteina-annotated-cif/scripts/check_annotated_cif.sh` on all four.

### Checker output, and why it is clean

The checker reports 1 error and 3 warnings per file. Running it on the repo's own
shipped `pmvk_3CH4__VHH_caplacizumab_cdr123.cif.zst` as a control produces the
**same** error and the same two substantive warnings, so none is introduced here:

| Finding | Control | Here |
| --- | --- | --- |
| `pipeline_error` AddChiralityFeature | yes | yes — Windows-only, see flag 10 |
| `component_ccd_atom_set_mismatch` (3 sentinels) | yes | yes — the CDR sentinels |
| `formal_charge_unknown` | yes (`+1: 28`) | yes (`+1: 24`, `-1: 17`) |
| `res_id_gap_check_skipped` | yes | **no** — index groups present, check ran |
| `struct_conn` records | 6 | **8** — the two added disulfides |

---

# Flags

## Resolved with you

1. **`A/54`, `A/56`, `A/96` do not exist.** Chain A of `target_crop.pdb` runs 151–290.
   Under a +150 offset (1-based index into the crop) they land on LEU 204, LEU 206,
   TYR 246 — residue names matching your LEU/LEU/TYR exactly, and all three appear in
   `CSS_epitope_residue_report.csv`. Confirmed: use 204/206/246. Note this makes the
   request *inconsistent with the template*, whose hotspots (`D166`, `Q241`, `Y246`)
   are written in auth numbering.

2. **Both `CD1` atoms are buried.** Shrake–Rupley per-atom SASA on the isolated target:

   | Residue | CD1 | CD2 |
   | --- | --- | --- |
   | LEU 204 | **0.00 Å²** | 26.92 Å² |
   | LEU 206 | **3.52 Å²** | 21.39 Å² |

   CD1/CD2 naming on leucine is largely arbitrary between structures. Switched to CD2.
   `hotspot_obedience` is an all-heavy-atom 4.5 Å nearest-neighbour criterion, and the
   template gates it at `min: 1.0` — all hotspots must be contacted — so a 0.00 Å² atom
   would have made the gate structurally unreachable for every design.

3. **No atom name was given for TYR 246.** You first chose CZ; its SASA is 0.50 Å²
   (buried). Switched to **CE1** (17.74 Å²), the exposed ring carbon, preserving the
   aromatic/hydrophobic intent. For reference OH is 43.97 Å² and was the template's atom.

4. **`target_crop.pdb` contains no binder** — one chain, 1134 atoms, target only. The
   template paired the target with the 1G8 nanobody as chain H. Resolved by using the
   four frameworks the repo ships in `examples/antibody_design/structures/vhh/`.

## Still open — read these before running

5. **The template CIF was never read.** `inputs:` in the template points at
   `/net/scratch/dlsee/projects/gates/scally/.../css_1g8pin_1G8_cdr123_L13-10-14_D166OD2_Q241NE2_Y246OH.cif`,
   which is not mounted on this machine. The conditioning above was reconstructed from
   the template YAML's prose description plus the repo's example inputs, **not copied
   from the actual file**. If that file's annotations differ from the examples', this
   build differs with them. Worth a diff on the cluster.

6. **"Same C_CRDGRP (backbone-only)" is ambiguous, and I chose the template's reading.**
   In the template and in the repo examples the *target* carries **absolute** backbone
   coordinates (group 0) and only the *framework* is internal (group 1). That is what
   was built. But `AGENTS.md` states the PPI default as *target* backbone-only
   **internal** coordinates, and `docs/how-to/input_prep/examples/condition_spec_ppi.yaml`
   does exactly that — target in `C_CRDGRP` group 1, pose free. The two differ in whether
   the target's global pose is pinned. Say the word and it is a one-line change in
   `build_annotated_cifs.py`.

7. **"Preserve bonds" had nothing to preserve, and that is itself the problem.**
   `target_crop.pdb` has no `CONECT`, `SSBOND` or `LINK` records. It does have four
   cysteines forming **two real disulfides** — CYS170–CYS195 and CYS226–CYS264, both at
   2.04 Å S–S. AtomWorks does **not** infer these from distance, so a literal
   "preserve what is there" build would have produced a target with two broken
   disulfides. They were added explicitly (`BondType.SINGLE`, verified < 2.5 Å) and
   appear as the 2 extra `_struct_conn` records. **This is a deliberate chemistry
   change relative to the input file** — reverse it if it was not what you wanted.

8. **The new epitope is markedly harder than the template's.** Measured on the final
   atom set, with the template's own numbers (from its YAML comments) alongside:

   | | This epitope (204 CD2 / 206 CD2 / 246 CE1) | Template (D166 OD2 / Q241 NE2 / Y246 OH) |
   | --- | --- | --- |
   | max angular spread from centroid | **57.7°** | 37.5° |
   | burial proxy (heavy atoms within 8 Å) | **78 / 84 / 67** | 34 / 40 / 48 |
   | max pairwise span | 13.77 Å | 19.87 Å |

   The span is comfortably inside a VHH paratope's ~20–25 Å reach, so the three atoms
   are co-reachable in principle. But they wrap 57.7° around the antigen rather than
   sitting on one flat face, and they are 2× more recessed than the template's triad.
   This is a groove, not a patch. Expect a lower `hotspot_adherence` pass rate than the
   D166/Q241/Y246 arm, and do not read a drop as a framework effect.

9. **CDR lengths are sampled here, not pinned — so this is not comparable to the
   length-pinned 1G8 run.** The template you named is the *length-pinned positive
   control* (`L13-10-14`). The repo example frameworks ship windows **13-13 / 9-10 /
   12-21**, which are neither 1G8's native 13/10/14 nor the lab template's authored
   11-16 / 8-12 / 8-22. One of the three (13-13) happens to be pinned; the other two
   sample. If you want a pinned arm, say which lengths and I will materialize them.

10. **A repo bug blocks `run_pipeline_validation` on Windows.**
    `src/rfproteina/transforms/chirality_feature.py:171` seeds with
    `np.random.randint(0, 2**32)`. The legacy NumPy API defaults to dtype `long`,
    which is int32 on Windows, so it raises `high is out of bounds for int32` and the
    pipeline check dies before validating anything. On Linux `long` is 64-bit and it
    passes, which is why the repo (classified POSIX-only) has not hit it. The repo's
    *own* example inputs fail the checker here for the same reason. `build_annotated_cifs.py`
    works around it with a scoped, documented shim so validation genuinely ran; the
    shim only seeds a shuffle of chirality frames and touches nothing written to disk.
    The real fix is one line (`np.random.default_rng()`), and is worth a PR.

11. **`examples/antibody_design/vhh.yaml` would throw your hotspots away.** It sets
    `generation.dataset.sampling_transforms: SampleNearbyHotspots`, which resamples
    hotspots per design. If you run these inputs through that config, the three atoms
    you chose are discarded. Per the example README, pass
    `generation.dataset.sampling_transforms=null` to retain them — or use a
    campaign-style config derived from the 1G8 template, which has no such transform.

12. **No generate config was produced.** Only CIFs were built. The template's `inputs`,
    `out_dir` and `run_name` are cluster paths, and its `metrics_tags`, checkpoint pins,
    `recycle_strategy: medium/front_loaded` and four pass criteria were left untouched.
    Its metric masks still work unchanged: `sa_mask('C_SEQ')` on chain H still separates
    framework from CDRs, and the hotspot count is still 3, so `hotspot_obedience >= 1.0`
    still means all three.

13. **Cosmetic: the binder sits in the PMVK pose.** Chain H keeps the coordinates it had
    against PMVK in the source example, so in a viewer the nanobody will not be docked
    onto the CSS target. Only intra-group geometry is conditioned and `C_CRD` is NaN
    across chain H, so the model samples placement from noise — exactly the situation
    the template YAML describes for its own framework.
