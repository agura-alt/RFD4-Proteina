# Design run history

Reconstructed 2026-11-19 from the output trees under `runs/`, the provenance headers in
`CSS/configs/*.yaml`, and git history. Updated 2026-11-23 with the Phase 3 results.

Provenance of the numbers below: Phase 1, Phase 2 and the native control were recomputed
locally from each run's `per_example_metrics.csv`. The Phase 3 arms have **not** been synced
to this machine — their figures are transcribed from `CSS/report_ipae.py` output run against
the output trees on Delta. Sync the CSVs before relying on anything beyond what is quoted here.

---

## Read this before reading any number

`oracle_ppi.min_ipae_coi` is AF3's `chain_pair_pae_min` passed through unscaled
(`src/rfproteina/metrics/oracles/oracle_metric.py`, `_COI_INTERFACE_METRICS`): the **minimum
predicted aligned error between the two chains, in Ångströms**, saturating near 31.5.

That makes the scale, not the ranking, the thing to read:

| Range | What AF3 is saying |
|---|---|
| 15–25 Å | No predicted interface. The binder's position relative to the antigen is unconstrained. |
| 5–12 Å | Sparse transition region; few designs land here. |
| ≤ 5 Å | A confidently predicted interface. |
| ≤ 1.2 Å | The campaign's `min_ipae` pass criterion. |

Two consequences that invalidate how earlier revisions of this file read the data:

1. **The bulk is not "poor designs."** It is the no-interface floor. A distribution centred at
   20 Å contains no interface predictions at all, so comparing arm means over that bulk
   compares how much nothing each arm produced.
2. **The 1.2 Å gate sits near the floor of AF3's reported range**, well inside the conventional
   ≤5 Å confident-interface threshold. A pass rate of 0/1000 against it is consistent with a
   campaign that produced confident interfaces at a useful rate. It did — see Phase 3.

The distributions are close to bimodal for this reason: a no-interface bulk, a sparse valley,
and a small confident mode. Read the ≤5 Å fraction, not the mean.

---

## Phase 1 — Exploratory VHH pilots (repo example targets)

Repo example antigens, not CSS. Establishes that the pipeline runs end to end.

| Run | Target / framework | n | min_ipae_coi best / mean | Gates |
|---|---|---|---|---|
| `runs/first` | first binder smoke | 2 | no oracle metrics | — |
| `runs/vhh_pilot` | mzb1 (7AAH) + pmvk (3CH4) × 4 VHH frameworks | 8 | 3.36 / 10.12 | 0 pass |
| `runs/vhh_pilot2` | same, re-run | 8 | 1.79 / 10.48 | 0 pass |

---

## Phase 2 — CSS campaign, epitope L204.CD2 / L206.CD2 / Y246.CE1

Four repo VHH frameworks (caplacizumab, ozoralizumab, vobarilizumab, h-NbBCII10) against
the CSS target crop. Inputs built by `CSS/build_annotated_cifs.py`; configs in `CSS/configs/`.

### Smoke tests, 2026-11-18

`runs/smoke_20261118_145742`, `runs/smoke_20261118_160816`, `runs/smoke_20261118_161549`,
`runs/submit_test`. Staged 1-design (stage0) → 12-design (stage1) runs on medium /
front-loaded recycling, 180 trunk passes. 12 designs at 15.8–24.5 Å — i.e. entirely within
the no-interface bulk, which at n=12 is the expected outcome and not evidence of anything.

### High-recycling arm — `runs/css_high`

`CSS/configs/css_cdr123_L204CD2_L206CD2_Y246CE1_high_generate.yaml`: high / front-loaded,
300 trunk passes, 100 designs (25 per framework).

| Framework | iPAE mean ± sd | best | hotspot | apo/holo | apo pTM |
|---|---|---|---|---|---|
| vobarilizumab | **18.27 ± 2.78** | 12.09 | **15/25** | 9/25 | 66.8 |
| ozoralizumab | 20.28 ± 3.61 | **8.96** | 9/25 | 3/25 | 54.0 |
| h-NbBCII10 | 20.56 ± 3.43 | 15.17 | 14/25 | 10/25 | 70.8 |
| caplacizumab | 21.32 ± 2.75 | 15.10 | 11/25 | 17/25 | 71.0 |

Totals: 0/100 on `min_ipae` and `ppi`; 49/100 hotspots, 39/100 apo/holo. **0/100 under 5 Å** —
no confident interface in this arm. At the ~2–3% rate Phase 3 later measured, 100 designs
expect ~2, so this is a depth limit, not a framework verdict.

Framework choice is the one manipulation that has moved anything: vobarilizumab separates
from caplacizumab at t ≈ 3.9, which is why it was carried forward alone into the depth probes.
The 8.96 design is an ozoralizumab artefact — apo pTM 37, so its binder barely folds alone.

### Native positive control — BROKEN, do not cite as calibration

`runs/native_control`, `_2`, `_3`, and `runs/native_control_annotated`
(`CSS/configs/css_native_control_evaluate.yaml`, input from
`CSS/build_native_control_annotated.py`). The crystallographic 1G8 complex scored through
the campaign's own metric definitions.

| Run | min_ipae_coi | iptm_coi | plddt_coi | ptm_coi | apo plddt_coi | apo ptm_coi |
|---|---|---|---|---|---|---|
| unannotated, plain `oracle` | 16.56 | 13 | — | — | — | — |
| annotated, `oracle_ppi` + `oracle_apo` | 14.19 | 17 | **52.1** | 51 | **39.6** | 31 |

A crystallographic complex should refold at iptm ~80+ and pLDDT ~90. pLDDT 52 on the chain of
interest, and 39.6 apo, means AF3 never folded the binder chain — so it was never in a position
to dock it. **This run is a broken control, not a calibrated one.** It is evidence that the
control input is wrong (chain mapping, annotation, or the target crop), and it says nothing
about whether the gate is reachable.

Earlier revisions of this file drew the opposite conclusion — that a known-good binder landing
inside the design distribution proved the gate miscalibrated. That reading is withdrawn. The
control never produced a folded binder, so it cannot bound what the gate can reach.

Rebuilding this control is the highest-value open item: without it, nothing anchors the metric.

---

## Phase 3 — Epitope and hotspot-sampling arms (on Delta, not synced here)

Design counts from `ls inference_samples | wc -l`; distributions from `CSS/report_ipae.py`.

| Arm | Epitope | Frameworks | Recycling | n | mean ± sd | min | ≤5 Å | ≤1.2 Å |
|---|---|---|---|---|---|---|---|---|
| `css_wallphilic_E165OE2_D166OD2_H239CE1_Q241NE2` | wall, 4 hotspots | 4 | medium/front_loaded | 600¹ | 20.12 ± 4.40 | 1.70 | **16 (2.7%)** | 0 |
| `css_6hot_E165OE2_D166OD2_L204CD2_L206CD2_Q241NE2_Y246CE1` | 6 hotspots | 4 | medium/front_loaded | 1000 | 20.20 ± 4.49 | 1.78 | **26 (2.6%)** | 0 |
| `css_vobarilizumab_1000_hotspotadd` | L204/L206/Y246 + 2–3 drawn | 1 (voba) | high/front_loaded | 1000 | 18.99 ± 4.17 | 1.27 | **17 (1.7%)** | 0 |

¹ Incomplete: 600 of 1000, i.e. six of ten array tasks. `skip_existing` is true, so resubmitting
to the same `out_dir` backfills only the gap.

**The headline: ~59 designs across these three arms predict a confident interface (≤5 Å).**
That is 1.7–2.7% per arm, in the neighbourhood of the ~1% published rate for hotspot-conditioned
binder design. The campaign produced leads; the 1.2 Å gate did not report them.

The mode is real and not a Gaussian tail. A normal(20.2, 4.49) at n=1000 has an expected minimum
of ~5.7; the observed minima are 1.27–1.78, and the histograms show a populated low cluster
separated from the bulk by a sparse valley.

### What did and did not move

- **Epitope: no effect.** wallphilic 20.12 ± 0.18 vs 6hot 20.20 ± 0.14 (SEs) — 0.08 apart, under
  half a standard error. Two different epitopes, one conditioning on six hotspots rather than
  four, giving superimposable distributions and the same ≤5 Å rate.
- **Hotspot resampling: no effect.** `hotspotadd` (18.99, n=1000) against the matched
  `css_high` vobarilizumab slice (18.27 ± 2.78, n=25 — same epitope, framework and recycling)
  is **+0.72 ± 0.57, t = 1.27**: not significant, and the sign is *worse*. n=25 is sufficient
  here because the baseline SE is 0.57.
- **`hotspotadd`'s lower mean is framework composition, not sampling.** It is vobarilizumab-only
  against four-framework mixtures, and vobarilizumab already ran ~1.2 Å below the panel mean.

### Comparing arms: what is and is not valid

- **`min`, `max` and the best-N tables are n-dependent** and must not be compared across arms of
  different size. `hotspotadd`'s 1.27 against wallphilic's 1.70 is substantially 1000 draws
  against 600. When wallphilic's last tasks land its minimum will likely fall to ~1.3–1.5 from
  the extra draws alone; that is not an improvement.
- **Means, medians and threshold fractions are safe** at these n (SE 0.14–0.18).
- Only **wallphilic vs 6hot** is a clean two-arm comparison — same frameworks, recycling and
  reps, with epitope the single variable. `hotspotadd` differs from both in epitope, framework
  count, recycling budget *and* sampling scheme, so its row in a three-way table is not
  interpretable as a hotspot-sampling result.
- Hotspot *pass rates* are not comparable across the sampling arms at all: `hotspot_obedience`
  is a fraction over 3, 5–6 and 2–3 conditioned atoms respectively.

### Not run

| Config | Status |
|---|---|
| `css_cdr123_..._hotspotadd / _hotspotswap` | **Never submitted.** Committed 2026-11-19 (`69882a8`); no output tree exists. 400 designs each if wanted. |
| `css_vobarilizumab_1000_generate` (fixed baseline) | Not on Delta. The matched 25-design comparison above stands in for the bulk question; a 1000-design version would only add tail resolution. |
| `css_vobarilizumab_1000_hotspotswap` | Not on Delta. Unconfirmed elsewhere — the listing that established this covered one cluster only. |
| `css_vobarilizumab_10000_generate` | Running off-Delta; may not complete. Partial output is still usable — every finished array task writes its own metrics, and `report_ipae.py` pools whatever is present. |
| `css_groove_E165OE2_L204CD2_L206CD2_Q241NE2` | Status unknown; not covered by the listings taken so far. |

---

## Through-line

The campaign produced ~59 designs with confidently predicted interfaces, at 1.7–2.7% per arm,
and none of them was visible through the `min_ipae ≤ 1.2` gate the campaign reported on.

Neither epitope choice nor hotspot resampling moved that rate. Framework choice is still the
only manipulation with a measurable effect. Depth is what surfaced the mode: it is invisible at
n=100 and clear at n=600–1000, which was the depth probes' purpose, and it is already answered
without the 10k.

Open, in priority order:

1. **Triage the ≤5 Å designs.** Needs the Phase 3 CSVs synced. Screen on `iptm_coi`,
   `plddt_coi`, `oracle_apo.apo.ptm_coi` (does the binder fold alone — this is what disqualified
   the 8.96 design), `has_clash`, hotspot obedience and CDR contact counts, then cluster with
   foldseek so distinct solutions are handed over rather than near-duplicates. Note that no
   design currently passes every gate: the 1.27 passes `ppi` and `antibody` with hotspot
   obedience 0.5; wallphilic's 1.70 passes hotspot and `antibody` but fails `ppi`.
2. **Rebuild the native control**, which is the only thing that would anchor the metric.
3. **Re-set the reporting threshold.** ≤5 Å is the defensible confident-interface cut; 1.2 Å
   as a *pass* criterion needs justification it has never had.
4. Backfill wallphilic's four missing array tasks (resubmit, same `out_dir`).

These are AF3-confidence leads. In-silico interface confidence is an optimistic signal and
routinely survives into designs that fail experimentally — worth ordering, not worth promising.
