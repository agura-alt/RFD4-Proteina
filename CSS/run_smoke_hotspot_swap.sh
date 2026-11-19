#!/usr/bin/env bash
# Smoke test for the HOTSPOT-SWAP arm of campaign css_cdr123_L204CD2_L206CD2_Y246CE1.
#
# 2-3 hotspots sampled from a patch around the three chosen atoms, INSTEAD of them.
# Its sibling CSS/run_smoke_hotspot_add.sh keeps the chosen three and adds the sampled
# atoms to them; plain CSS/run_smoke.sh is the fixed three.
#
# Same two stages, same inputs, same gates as CSS/run_smoke.sh. The one difference is
# the config: css_cdr123_L204CD2_L206CD2_Y246CE1_hotspotswap_generate.yaml adds
# `generation.dataset.sampling_transforms` with a bare SampleNearbyHotspots, which
# REPLACES the C_HOT seeds and excludes the seed residues (hotspot_sampling.py:87).
#
# READ THIS BEFORE USING THE RESULTS. A/204 LEU, A/206 LEU and A/246 TYR are NOT
# conditioned on any design in this arm; they only decide where the patch sits. The
# designs this produces are not hits on the selected epitope, and nothing in the
# configured metrics measures contact with it. Two further things the numbers will
# not tell you:
#   * the conditioned count is 2 or 3, and hotspot_adherence is a FRACTION gated at
#     >= 1.0, so it is EASIER here than the fixed arm's 3-of-3 and much easier than
#     the add arm's 5-of-5. A higher pass rate may be nothing but the 2-hotspot draws.
#   * every design conditions on a different surrounding patch, so a failure is
#     ambiguous between "this framework cannot dock" and "that draw was a bad patch".
#
# Run from the repository root, with the rfproteina environment active and network
# access to the Hugging Face Hub.

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
export PROJECT_ROOT="$PWD"

CONFIG_NAME=css_cdr123_L204CD2_L206CD2_Y246CE1_hotspotswap_generate

echo "############ Preflight: hotspot patch fit, no GPU ############"
# Does the resampling actually fit on these inputs, and did the sampling_transforms
# block compose at all? A mis-addressed key is silently inert: generate.py:204 finds
# no transform, the run falls back to the fixed seed hotspots -- which in THIS arm
# means it quietly becomes the baseline arm -- and every metric afterwards looks
# healthy. This is the only step that catches that.
python CSS/preflight_hotspot_patch.py 25 --config "$CONFIG_NAME"

# Everything else -- the metric preflight, stage 0, stage 1, the design-count guards
# and the column audit -- is CSS/run_smoke.sh, driven onto this config. Recycling is
# unchanged from the medium arm (PPI medium/front_loaded, 180 passes), so the
# preflight's pass expectation stays at its default.
echo
SMOKE_CONFIG_NAME="$CONFIG_NAME" \
SMOKE_EXPECTED_PASSES=180 \
SMOKE_LABEL=hotspotswap \
  bash CSS/run_smoke.sh
