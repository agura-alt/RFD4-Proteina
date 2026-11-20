# Design run history

Reconstructed 2026-11-19 from the output trees under `runs/`, the provenance headers in
`CSS/configs/*.yaml`, and git history. Numbers below were recomputed from each run's
`per_example_metrics.csv`, not copied from prose.

The campaign gate is `oracle_ppi.min_ipae_coi <= 1.2` (pass criterion `min_ipae`).
`CSS/report_ipae.py` prints the distribution behind that pass/fail, which matters here:
every failure so far has been an order of magnitude out, not a near miss.

---

## Phase 1 — Exploratory VHH pilots (repo example targets)

Repo example antigens, not CSS. Establishes that the pipeline runs end to end.

| Run | Target / framework | n | min_ipae_coi best / mean | Gates |
|---|---|---|---|---|
| `runs/first` | first binder smoke | 2 | no oracle metrics | — |
| `runs/vhh_pilot` | mzb1 (7AAH) + pmvk (3CH4) × 4 VHH frameworks | 8 | 3.36 / 10.12 | 0 pass |
| `runs/vhh_pilot2` | same, re-run | 8 | 1.79 / 10.48 | 0 pass |

The best iPAE in the entire history (1.79) comes from here — against an easy target, and
still short of the 1.2 gate.

---

## Phase 2 — CSS campaign, epitope L204.CD2 / L206.CD2 / Y246.CE1

Four repo VHH frameworks (caplacizumab, ozoralizumab, vobarilizumab, h-NbBCII10) against
the CSS target crop. Inputs built by `CSS/build_annotated_cifs.py`; configs in
`CSS/configs/`.

### Smoke tests, 2026-11-18

`runs/smoke_20261118_145742`, `runs/smoke_20261118_160816`, `runs/smoke_20261118_161549`,
`runs/submit_test`. Staged 1-design (stage0) → 12-design (stage1) runs on medium /
front-loaded recycling, 180 trunk passes.

- 12 designs at iPAE 15.8–24.5, **0/12** on the 1.2 gate, 7/12 hotspots, 4/12 apo/holo.
- `submit_test` (4 designs, one per framework) exercised the SLURM submission path;
  one design reached 6.23.

### High-recycling arm — `runs/css_high`

`CSS/configs/css_cdr123_L204CD2_L206CD2_Y246CE1_high_generate.yaml`: high / front-loaded,
300 trunk passes, 100 designs (25 per framework). The campaign's main dataset.

| Framework | iPAE mean ± sd | best | hotspot | apo/holo | apo pTM |
|---|---|---|---|---|---|
| vobarilizumab | **18.27 ± 2.72** | 12.09 | **15/25** | 9/25 | 66.8 |
| ozoralizumab | 20.28 ± 3.61 | **8.96** | 9/25 | 3/25 | 54.0 |
| h-NbBCII10 | 20.56 ± 3.43 | 15.17 | 14/25 | 10/25 | 70.8 |
| caplacizumab | 21.32 ± 2.75 | 15.10 | 11/25 | 17/25 | 71.0 |

Totals: **0/100** on `min_ipae` and `ppi`; 49/100 hotspots, 39/100 apo/holo.

Two readings:

- Tripling the recycling budget (180 → 300 passes) did not move the gate.
- The 8.96 best is an ozoralizumab design that fails hotspots with apo pTM 37 — its binder
  barely folds alone, so it reads as an artefact rather than a lead. Ozoralizumab was
  0/25 on the joint gates with the worst apo folding of the panel. Vobarilizumab is the
  only framework whose interface error separates from the others beyond noise
  (t ≈ 3.9 against caplacizumab), which is why it was carried forward alone.

### Native positive control

`runs/native_control`, `_2`, `_3`, and `runs/native_control_annotated`
(`CSS/configs/css_native_control_evaluate.yaml`, input from
`CSS/build_native_control_annotated.py`). The crystallographic 1G8 complex scored through
the campaign's own metric definitions.

| Run | min_ipae_coi | iptm | hotspot | apo/holo |
|---|---|---|---|---|
| unannotated, plain `oracle` | 16.56 | 0.13 | — | — |
| annotated, `oracle_ppi` + `oracle_apo` | **14.19** | 17 | pass | fail |

This is the most important result in the file. A binder known to work on this antigen lands
*inside* the design distribution (designs: 15.8–24.5) and misses the 1.2 gate by more than
an order of magnitude. Until that is resolved, a failing design and an unreachable gate are
indistinguishable — the gate, not the designs, is what is miscalibrated.

---

## Phase 3 — Configured, no outputs on this machine

Committed 2026-11-19 but with no results synced back here. Treat as not yet run unless the
cluster says otherwise.

| Config | Designs | Notes |
|---|---|---|
| `css_vobarilizumab_1000_generate.yaml` | 1,000 | vobarilizumab only, 300 passes, `gen_njobs` 10, compressed outputs |
| `css_vobarilizumab_1000_hotspotadd_generate.yaml` | 1,000 | 3 chosen atoms + 2–3 from a patch, redrawn per replicate |
| `css_vobarilizumab_1000_hotspotswap_generate.yaml` | 1,000 | 2–3 patch atoms *instead of* the chosen three |
| `css_vobarilizumab_10000_generate.yaml` | 10,000 | ~29 s/design ≈ 80 GPU-h, ~22 GB pre-compression, `gen_njobs` 200 |
| `css_cdr123_..._hotspotadd / _hotspotswap_generate.yaml` | 100 each | the same two hotspot arms on the four-framework campaign |
| `css_groove_E165OE2_L204CD2_L206CD2_Q241NE2_generate.yaml` | 1,000 | four-hotspot groove-spanning epitope, all four frameworks at 250 reps, medium/front_loaded |

Compare the hotspot arms on **iPAE, not on the hotspot gate** — `hotspot_obedience` is a
fraction over 3, 5–6 and 2–3 conditioned atoms respectively, so the denominators differ.

What the 10,000 probe is for: detection. 100 designs resolve ~5%, 1,000 resolve ~0.5%,
10,000 resolve ~0.05%. The published rate for hotspot-conditioned binder design is ~1%, so
10,000 either produces ~100 passers or excludes that rate outright. 112 designs so far have
produced nothing below 8.96, and at that size a 1% mode contributes about one design and is
indistinguishable from the tail of a unimodal bulk.

---

## Through-line

Across ~140 locally recorded CSS designs, `min_ipae <= 1.2` has never been passed, and the
native control says it may not be passable as written. Hotspot adherence is the one gate
that moves with the framework choice. The open question the depth probes exist to answer is
whether a passing subpopulation exists at all.
