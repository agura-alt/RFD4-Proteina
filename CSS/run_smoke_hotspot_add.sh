#!/usr/bin/env bash
# Smoke test for the HOTSPOT-ADD arm of campaign css_cdr123_L204CD2_L206CD2_Y246CE1.
#
# The three chosen hotspots PLUS 2-3 sampled from around them, per replicate.
# Its sibling CSS/run_smoke_hotspot_swap.sh drops the chosen three and conditions on
# the sampled atoms alone; plain CSS/run_smoke.sh is the fixed three.
#
# Same two stages, same inputs, same gates as CSS/run_smoke.sh. The one difference is
# the config: css_cdr123_L204CD2_L206CD2_Y246CE1_hotspotadd_generate.yaml adds
# `generation.dataset.sampling_transforms` with SampleNearbyHotspots wrapped in
# KeepSeedHotspots, so every design is conditioned on the three selected atoms PLUS
# 2-3 more drawn from a patch around them, redrawn per replicate.
#
# The wrapper is there because SampleNearbyHotspots on its own REPLACES the seeds and
# excludes the seed residues (hotspot_sampling.py:87) -- A/204 LEU, A/206 LEU and
# A/246 TYR would set where the patch sits and then be dropped. KeepSeedHotspots
# (src/rfproteina/transforms/sampling/css_seeded_hotspots.py) ORs them back in.
# NOTE that module's directory is git-ignored (.gitignore:9); the module itself was
# force-added so it travels with the config. Anything further added beside it needs the
# same treatment.
#
# READ THIS BEFORE USING THE RESULTS. Two things the numbers will not tell you:
#   * the conditioned count is 5 or 6, and hotspot_adherence is a FRACTION gated at
#     >= 1.0, so it now demands ALL of them within 4.5 A rather than all 3. That is a
#     strictly harder gate than the fixed arm's, and designs that satisfy the chosen
#     epitope perfectly well will fail it. Score the three seeds separately offline.
#   * every design conditions on a different surrounding patch, so a failure is
#     ambiguous between "this framework cannot dock" and "that draw was a bad patch".
#
# Run from the repository root, with the rfproteina environment active and network
# access to the Hugging Face Hub.

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
export PROJECT_ROOT="$PWD"

echo "############ Preflight: hotspot patch fit, no GPU ############"
# Does the resampling actually fit on these inputs, and did the sampling_transforms
# block compose at all? A mis-addressed key is silently inert: generate.py:204 finds
# no transform, the run falls back to the fixed seed hotspots, and every metric
# afterwards looks healthy. This is the only step that catches that.
python CSS/preflight_hotspot_patch.py 25 --config css_cdr123_L204CD2_L206CD2_Y246CE1_hotspotadd_generate

# Everything else -- the metric preflight, stage 0, stage 1, the design-count guards
# and the column audit -- is CSS/run_smoke.sh, driven onto this config. Recycling is
# unchanged from the medium arm (PPI medium/front_loaded, 180 passes), so the
# preflight's pass expectation stays at its default.
echo
SMOKE_CONFIG_NAME=css_cdr123_L204CD2_L206CD2_Y246CE1_hotspotadd_generate \
SMOKE_EXPECTED_PASSES=180 \
SMOKE_LABEL=hotspotadd \
  bash CSS/run_smoke.sh
