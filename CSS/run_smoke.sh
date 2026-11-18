#!/usr/bin/env bash
# Smoke test for campaign css_cdr123_L204CD2_L206CD2_Y246CE1.
#
# Two stages, deliberately. They fail for different reasons and cost very
# different amounts, so running them separately stops you paying for AF3 to
# discover a checkpoint typo.
#
#   Stage 0  1 design, 1 input      -- does the run start at all?
#   Stage 1  12 designs, 4 inputs   -- does every metric column appear?
#
# Neither stage says anything about design quality. 12 designs cannot estimate a
# pass rate; they exist to prove the gates are determinate before you spend 400.
#
# Run from the repository root, with the rfproteina environment active and
# DATA_PATH pointing at the checkpoints.

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
export PROJECT_ROOT="$PWD"

CONFIG_DIR="$PROJECT_ROOT/CSS"
CONFIG_NAME="css_cdr123_L204CD2_L206CD2_Y246CE1_generate"
INPUT_DIR="$PROJECT_ROOT/CSS/inputs/cifs_L204CD2_L206CD2_Y246CE1"
STAMP="$(date +%Y%m%d_%H%M%S)"

# Fresh output directories every time. `skip_existing` defaults to true, so a
# reused directory silently keeps the previous run's structures and you end up
# auditing stale results (AGENTS.md: use fresh output directories when changing
# settings).
STAGE0_DIR="$PROJECT_ROOT/runs/smoke_${STAMP}/stage0"
STAGE1_DIR="$PROJECT_ROOT/runs/smoke_${STAMP}/stage1"

# Shared overrides.
#   seed            pinned so a failure is reproducible while you debug it
#   batch_size=1    AGENTS.md: test the longest input at batch size 1 first
#   verbose         the diagnostics that say WHY a metric produced nothing
#   prevalidate     fail fast on a bad input rather than mid-batch
#   save_chain_of_interest
#                   needed for foldseek clustering later, and cannot be
#                   recovered without re-running; turn it on now so the flag is
#                   exercised and the production run is not a surprise
COMMON=(
  --config-path "$CONFIG_DIR"
  --config-name "$CONFIG_NAME"
  batch_size=1
  seed=42
  verbose=true
  prevalidate_inputs=true
  save_chain_of_interest=true
)

# `set -e` aborts on a non-zero exit, which gates the stages below. It does NOT
# catch this campaign's actual signature failure: exiting 0 having produced
# nothing. A mistyped input path, or an input glob that matches no file, yields a
# clean exit and zero structures -- and without this guard the script would walk
# straight into stage 1 and spend 24 AF3 calls on it.
assert_designs() {
    local dir="$1" min="$2" label="$3" n
    n=$(find "$dir" -type f \( -name '*.cif' -o -name '*.cif.zst' \) 2>/dev/null | wc -l)
    if [ "$n" -lt "$min" ]; then
        echo "FAIL: $label produced $n structure file(s) under $dir, expected at least $min" >&2
        echo "      The job exited 0, so this would otherwise have looked like success." >&2
        return 1
    fi
    echo "ok: $label produced $n structure file(s)"
}

echo "############ Preflight: config + metric selection, no GPU ############"
# Static check: composes the config and resolves which metrics the tag list
# selects, using the same rule the runtime applies. Costs seconds and catches a
# mis-addressed recycle_strategy or a dead tag before anything is scheduled.
python CSS/preflight_metric_selection.py

echo
echo "############ Stage 0: plumbing, 1 design ############"
# Largest of the four inputs (1781 atoms; the others are within 11 atoms, so this
# is the worst case for memory). The 2026-11-04 failures died ~36 s in on a
# checkpoint error -- that class of bug surfaces here, before any AF3 call.
rfd4-generate "${COMMON[@]}" \
  inputs="$INPUT_DIR/css_VHH_caplacizumab_cdr123_L204CD2_L206CD2_Y246CE1.cif" \
  num_replicates=1 \
  out_dir="$STAGE0_DIR"

# Gate stage 1 on stage 0 having really made something, not just on exit 0.
assert_designs "$STAGE0_DIR" 1 "stage 0"

echo
echo "############ Stage 1: metric coverage, 12 designs ############"
# num_replicates is PER INPUT, so 3 x 4 frameworks = 12 designs, matching the
# size of this campaign's own 2026-11-04 smoke run. Three per input rather than
# one so that a single degenerate design cannot make a live column look absent.
rfd4-generate "${COMMON[@]}" \
  inputs="$INPUT_DIR" \
  num_replicates=3 \
  out_dir="$STAGE1_DIR"

# 12 designs expected (3 replicates x 4 frameworks). Fewer means inputs were
# skipped, which the audit below would not notice -- it reads whatever rows exist.
assert_designs "$STAGE1_DIR" 12 "stage 1"

echo
echo "############ Audit ############"
# Exit status is not the test. The 2026-11-04 run exited 0 with every hotspot
# gate indeterminate.
python CSS/audit_smoke_metrics.py "$STAGE1_DIR"
