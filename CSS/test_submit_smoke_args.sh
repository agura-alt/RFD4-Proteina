#!/usr/bin/env bash
# Unit test for the SMOKE argument logic in CSS/submit_generate.sh.
#
# That logic decides how many designs a batch job makes, so getting it wrong either wastes an
# allocation on 400 designs that were meant to be 12, or silently keeps a user's own override
# from taking effect. None of it needs a GPU, a checkpoint or Slurm, so it is tested here.
#
# The functions are extracted from the submit script between its `# >>> smoke args` and
# `# <<< smoke args` markers and sourced, so this tests the real code rather than a copy.
#
# Run from the repository root:  bash CSS/test_submit_smoke_args.sh

set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
SCRIPT="CSS/submit_generate.sh"
BLOCK="$(sed -n '/^# >>> smoke args/,/^# <<< smoke args/p' "$SCRIPT")"
if [ -z "$BLOCK" ]; then
    echo "FAIL: could not find the smoke-args markers in $SCRIPT" >&2
    exit 1
fi
GUARD="$(sed -n '/^# >>> guards/,/^# <<< guards/p' "$SCRIPT")"
if [ -z "$GUARD" ]; then
    echo "FAIL: could not find the guard markers in $SCRIPT" >&2
    exit 1
fi
# Both blocks are pure function definitions. CONFIG_NAME and CONFIG_FILE are the only outside
# names the smoke block reads; the guard reads NJOBS and the two SLURM_ARRAY_* variables.
CONFIG_NAME="css_test_arm"
PWD_SAVED="$PWD"
eval "$BLOCK"
eval "$GUARD"

# Fixtures for count_inputs: a four-input arm (a directory of CIFs) and a single-input
# depth probe. Real files, because count_inputs resolves the path on disk.
FIXTURES="$(mktemp -d)"
trap 'rm -rf "$FIXTURES"' EXIT
mkdir -p "$FIXTURES/cifs"
for i in 1 2 3 4; do : > "$FIXTURES/cifs/framework_$i.cif"; done
: > "$FIXTURES/cifs/notes.txt"  # must not be counted
printf 'inputs: %s/cifs\n' "$FIXTURES" > "$FIXTURES/four.yaml"
printf 'inputs: %s/cifs/framework_1.cif\n' "$FIXTURES" > "$FIXTURES/one.yaml"
printf 'inputs: ${oc.env:PROJECT_ROOT}/cifs\n' > "$FIXTURES/interp.yaml"

failures=0
check() {
    local label="$1" expected="$2" actual="$3"
    if [ "$expected" = "$actual" ]; then
        printf 'ok   %-56s %s\n' "$label" "$actual"
    else
        printf 'FAIL %-56s expected [%s] got [%s]\n' "$label" "$expected" "$actual"
        failures=$((failures + 1))
    fi
}

# --- smoke_enabled reads SMOKE ------------------------------------------------
for value in 1 true TRUE yes YES; do
    SMOKE="$value" && smoke_enabled && r=on || r=off
    check "smoke_enabled with SMOKE=$value" "on" "$r"
done
for value in 0 "" no false off bogus; do
    SMOKE="$value" && { smoke_enabled && r=on || r=off; }
    check "smoke_enabled with SMOKE='$value'" "off" "$r"
done
unset SMOKE
smoke_enabled && r=on || r=off
check "smoke_enabled with SMOKE unset" "off" "$r"

# --- does --array agree with NJOBS? -------------------------------------------
# Each case is a submission that Slurm accepts and that exits 0 while generating the wrong
# set of designs, so the guard's exit status is the whole point. "" means the variable is
# unset, which is what a non-array submission looks like.
guard_rc() {
    (
        export NJOBS="$1"
        if [ -n "$2" ]; then export SLURM_ARRAY_TASK_COUNT="$2"; else unset SLURM_ARRAY_TASK_COUNT; fi
        if [ -n "$3" ]; then export SLURM_ARRAY_TASK_ID="$3"; else unset SLURM_ARRAY_TASK_ID; fi
        check_array_njobs >/dev/null 2>&1
        echo $?
    )
}
check "array 10, NJOBS 10: accepted"            "0" "$(guard_rc 10 10 0)"
check "array 10, NJOBS 1: refused"              "1" "$(guard_rc 1 10 0)"
check "array 4, NJOBS 10: refused"              "1" "$(guard_rc 10 4 0)"
check "NJOBS 10, no array at all: refused"      "1" "$(guard_rc 10 "" "")"
check "NJOBS 1, no array: accepted"             "0" "$(guard_rc 1 "" "")"
check "NJOBS 1, single-task array: accepted"    "0" "$(guard_rc 1 1 0)"
check "NJOBS unset defaults to 1: accepted"     "0" "$(NJOBS=1 guard_rc 1 "" "")"
# The message has to name the fix, since the log is all the user gets back from a dead job.
msg="$( (export NJOBS=10; unset SLURM_ARRAY_TASK_ID SLURM_ARRAY_TASK_COUNT; check_array_njobs 2>&1 >/dev/null) )"
case "$msg" in
    *"--array=0-9"*) check "no-array message names the fix" "yes" "yes" ;;
    *) check "no-array message names the fix" "yes" "no: $msg" ;;
esac

# --- was a script setting passed as a Hydra override? -------------------------
# The 10,000-design probe died this way: CONFIG_NAME after the script name reaches Hydra,
# which has no such key. A real Hydra override with a similar-looking name must still pass.
args_rc() { check_env_in_args "$@" >/dev/null 2>&1; echo $?; }
check "CONFIG_NAME as an override: refused"   "1" "$(args_rc CONFIG_NAME=css_vobarilizumab_10000_generate)"
check "NJOBS as an override: refused"         "1" "$(args_rc NJOBS=10)"
check "SMOKE as an override: refused"         "1" "$(args_rc SMOKE=1)"
check "REPO as an override: refused"          "1" "$(args_rc REPO=/scratch/x)"
check "refused even when not first"           "1" "$(args_rc num_replicates=4 CONFIG_NAME=x)"
check "real hydra overrides: accepted"        "0" "$(args_rc num_replicates=4 out_dir=/scratch/x seed=7)"
check "no arguments at all: accepted"         "0" "$(args_rc)"
check "config_name lowercase is hydra's"      "0" "$(args_rc config_name=x)"
check "a key merely containing it: accepted"  "0" "$(args_rc my_CONFIG_NAME=x)"
# The message must point at the env-prefix form, not at Hydra's '+KEY=' suggestion, which
# would append the key and silently run the default arm.
msg="$(check_env_in_args CONFIG_NAME=css_vobarilizumab_10000_generate 2>&1 >/dev/null)"
case "$msg" in
    *"BEFORE sbatch"*) check "message points at the env prefix" "yes" "yes" ;;
    *) check "message points at the env prefix" "yes" "no: $msg" ;;
esac

# --- how many inputs does the config name? ------------------------------------
CONFIG_FILE="$FIXTURES/four.yaml"
check "four-input arm: counted" "4" "$(count_inputs)"
CONFIG_FILE="$FIXTURES/one.yaml"
check "single-input depth probe: counted" "1" "$(count_inputs)"
CONFIG_FILE="$FIXTURES/interp.yaml"
check "PROJECT_ROOT is interpolated" "4" "$(PROJECT_ROOT="$FIXTURES" count_inputs)"
CONFIG_FILE="$FIXTURES/missing.yaml"
check "unreadable config: unknown" "0" "$(count_inputs)"
CONFIG_FILE=""
check "no config file: unknown" "0" "$(count_inputs)"

# The expectation follows the count, so a single-input arm is not failed for producing 3.
CONFIG_FILE="$FIXTURES/one.yaml"
build_smoke_args
check "single-input arm: expectation is 3" "3" "$SMOKE_EXPECTED"
CONFIG_FILE=""
build_smoke_args
check "unknown input count: expect only 'something'" "1" "$SMOKE_EXPECTED"

# --- defaults when the caller passes nothing ----------------------------------
CONFIG_FILE="$FIXTURES/four.yaml"
build_smoke_args
check "no args: replicates" "num_replicates=3" "$(printf '%s\n' "${SMOKE_ARGS[@]}" | grep '^num_replicates=')"
check "no args: batch size" "batch_size=1" "$(printf '%s\n' "${SMOKE_ARGS[@]}" | grep '^batch_size=')"
check "no args: seed pinned" "seed=42" "$(printf '%s\n' "${SMOKE_ARGS[@]}" | grep '^seed=')"
check "no args: verbose" "verbose=true" "$(printf '%s\n' "${SMOKE_ARGS[@]}" | grep '^verbose=')"
check "no args: prevalidate" "prevalidate_inputs=true" "$(printf '%s\n' "${SMOKE_ARGS[@]}" | grep '^prevalidate_inputs=')"
check "no args: expected designs" "12" "$SMOKE_EXPECTED"
check "no args: out_dir is fresh and named" \
    "$PWD_SAVED/runs/smoke_batch_css_test_arm_local_0" "$SMOKE_OUT_DIR"
check "no args: out_dir is passed to hydra" \
    "out_dir=$PWD_SAVED/runs/smoke_batch_css_test_arm_local_0" \
    "$(printf '%s\n' "${SMOKE_ARGS[@]}" | grep '^out_dir=')"

# --- the caller's own overrides win -------------------------------------------
build_smoke_args num_replicates=6
check "explicit replicates: not re-added" "" "$(printf '%s\n' "${SMOKE_ARGS[@]}" | grep '^num_replicates=')"
check "explicit replicates: no design expectation" "0" "$SMOKE_EXPECTED"
check "explicit replicates: other defaults survive" "batch_size=1" \
    "$(printf '%s\n' "${SMOKE_ARGS[@]}" | grep '^batch_size=')"

build_smoke_args out_dir=/scratch/mine
check "explicit out_dir: not re-added" "" "$(printf '%s\n' "${SMOKE_ARGS[@]}" | grep '^out_dir=')"
check "explicit out_dir: audited path follows it" "/scratch/mine" "$SMOKE_OUT_DIR"

build_smoke_args batch_size=4 seed=7 verbose=false
for key in batch_size seed verbose; do
    check "explicit $key: not re-added" "" "$(printf '%s\n' "${SMOKE_ARGS[@]}" | grep "^$key=")"
done

# A value containing '=' must survive, and a key that merely starts the same must not match.
build_smoke_args 'out_dir=/scratch/a=b' num_replicates_extra=9
check "out_dir value with '=' kept whole" "/scratch/a=b" "$SMOKE_OUT_DIR"
check "prefix-only key does not shadow" "num_replicates=3" \
    "$(printf '%s\n' "${SMOKE_ARGS[@]}" | grep '^num_replicates=')"

# --- sharding lowers the design expectation -----------------------------------
# Under --array the 12 designs are split across tasks; expecting 12 per task would fail
# every one of them.
NJOBS=4 build_smoke_args
check "sharded: expectation drops to 'something'" "1" "$SMOKE_EXPECTED"
NJOBS=1 build_smoke_args
check "unsharded: full expectation" "12" "$SMOKE_EXPECTED"

# --- Slurm ids reach the directory name ---------------------------------------
SLURM_JOB_ID=9988 SLURM_ARRAY_TASK_ID=2 build_smoke_args
check "job and array id in out_dir" \
    "$PWD_SAVED/runs/smoke_batch_css_test_arm_9988_2" "$SMOKE_OUT_DIR"

echo
if [ "$failures" -eq 0 ]; then
    echo "all checks passed"
    exit 0
fi
echo "$failures check(s) failed"
exit 1
