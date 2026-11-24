#!/bin/bash
# Submit a CSS generation run on Delta WITHOUT scripts/launch_generate.sh.
#
# Why not the launcher
# --------------------
# The launcher pre-spawns an AF3 daemon and points the oracle at its Unix socket via
# AF3_LAUNCHER_SOCKET_TEMPLATE. On Delta (2026-11-18, job 24147097) the daemon's srun step
# never produced a socket -- af3_daemons.log was empty -- and the run died at
# `sock.connect(daemon_sock)` with FileNotFoundError, after generating all four designs.
# The likely cause is that the daemon is an overlapping srun step requesting no GPU, which a
# GPU-type account may not submit ("Jobs requesting no GPU resources must be submitted under
# CPU type accounts" -- the same rule that rejects the launcher's concat job here). It also
# runs under /usr/bin/python3 rather than the venv.
#
# `_dispatch_fold` (af3_oracle.py) uses the daemon only when AF3_LAUNCHER_SOCKET_TEMPLATE is
# set, and otherwise calls apptainer directly. That direct path is what interactive runs on
# this cluster have been using all along. There is no launcher flag to turn the daemon off --
# the block is emitted unconditionally at launch_generate.sh:452 -- so this script submits
# the job itself and simply never sets the variable.
#
# What is given up: the launcher's job-array fan-out and its dependent aggregation job. The
# array is available here via --array (see below); aggregation is a separate command, and the
# launcher's concat job was already being rejected by this allocation anyway:
#   python src/rfproteina/inference_utils/aggregate_metrics.py \
#       --input <out_dir>/inference_metrics --cluster_with_foldseek
#
# Usage, from the repository root. Any extra arguments are passed through to rfd4-generate as
# Hydra overrides:
#
#   sbatch CSS/submit_generate.sh                        # the config's own settings
#   sbatch CSS/submit_generate.sh num_replicates=1 out_dir=$PWD/runs/submit_test
#   sbatch --time=12:00:00 CSS/submit_generate.sh        # override a directive
#
# Every task also writes its wall clock to slurm_run_outputs/timings/<config>_<jobid>.json
# (one JSON object per array task, written even when the task fails). To total an arm:
#   cat slurm_run_outputs/timings/<config>_*.json > timings.jsonl
# Slurm's own view of the same thing, which this does not replace:
#   sacct -j <jobid> -X --format=JobID,State,Submit,Start,End,Elapsed,ElapsedRaw
#
# Which arm: CONFIG_NAME picks the config in CSS/configs. sbatch exports the submitting
# environment by default (--export=ALL), so prefixing the variable is enough. The arms write
# to different out_dirs already, so they do not collide; give each its own job name and log
# files or three runs land in one %A-numbered pile:
#
#   C=css_cdr123_L204CD2_L206CD2_Y246CE1
#   for arm in _generate _hotspotadd_generate _hotspotswap_generate; do
#       CONFIG_NAME="$C$arm" sbatch --job-name="css$arm" \
#           --output="slurm_run_outputs/css${arm}_%A_%a.out" \
#           --error="slurm_run_outputs/css${arm}_%A_%a.err" \
#           CSS/submit_generate.sh
#   done
#
# The hotspot arms draw their conditioning per replicate, so `seed` decides which atoms each
# design gets. The configs do not pin one; pass seed=<n> if the draws must be reproducible.
# SMOKE=1 pins seed=42 for you.
#
# The 1000-design depth probes are the same three arms on ONE framework at the high
# recycling budget, and they take --array because gen_njobs is 10:
#
#   V=css_vobarilizumab_1000
#   for arm in _generate _hotspotadd_generate _hotspotswap_generate; do
#       NJOBS=10 CONFIG_NAME="$V$arm" sbatch --array=0-9 --job-name="voba$arm" \
#           --output="slurm_run_outputs/voba${arm}_%A_%a.out" \
#           --error="slurm_run_outputs/voba${arm}_%A_%a.err" \
#           CSS/submit_generate.sh
#   done
#
# Size: a plain submission runs the arm as configured -- 100 replicates x 4 inputs = 400
# designs and 800 AF3 refolds for the four-framework arms, or 1000 designs and 2000 refolds
# for a depth probe. SMOKE=1 runs stage 1 of CSS/run_smoke.sh instead (3 replicates PER
# INPUT -- 12 designs on a four-framework arm, 3 on a single-input probe -- batch size 1,
# verbose, prevalidated, a fresh out_dir), then gates on the design count and runs
# CSS/audit_smoke_metrics.py on the result:
#
#   SMOKE=1 CONFIG_NAME="${C}_hotspotadd_generate" sbatch --time=01:00:00 CSS/submit_generate.sh
#
# Shorten --time to match: the 4h default is sized for the production run. Any override you
# pass yourself wins over the smoke defaults, so `SMOKE=1 ... num_replicates=6` gives six.
#
# Resuming: submit again with the SAME out_dir. `skip_existing` defaults to true, so designs
# already written are skipped, and refolds whose confidence JSON exists are not repeated. That
# makes a short --time plus a resubmit a reasonable alternative to one long reservation. Note
# this is resuming an identical run; reusing an out_dir after CHANGING settings is the thing
# run_smoke.sh warns about, because the audit would then read stale structures.
#
# Fanning out: add --array=0-3 and set NJOBS=4 below so job_id/gen_njobs match the array size.

#SBATCH --job-name=css-generate
#SBATCH --account=bjua-delta-gpu
#SBATCH --partition=gpuA100x4
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --time=04:00:00
#SBATCH --signal=SIGUSR1@600
#SBATCH --output=slurm_run_outputs/css_generate_%A_%a.out
#SBATCH --error=slurm_run_outputs/css_generate_%A_%a.err

set -euo pipefail

# Wall-clock timing. Captured here, before any setup work, so the recorded elapsed covers the
# whole task -- venv activation, preflight, generation, the smoke gate -- and not just the
# model call. The matching report is written by timing_report() below, from an EXIT trap, so a
# task that fails or is killed still records what it spent. CLUSTER_LABEL only names the
# cluster in that report; it is deliberately NOT called CLUSTER, which generate.py's
# load_dotenv() reads to pick .env.<name> and which this script unsets further down.
TIMING_T0="$(date +%s)"
TIMING_START_ISO="$(date -Is)"
CLUSTER_LABEL="${CLUSTER_LABEL:-delta}"

REPO="${REPO:-/u/agurajapu/RFD4-Proteina}"
CONFIG_NAME="${CONFIG_NAME:-css_cdr123_L204CD2_L206CD2_Y246CE1_high_generate}"
NJOBS="${NJOBS:-1}"   # must equal the --array size when fanning out

# >>> guards -- extracted and unit-tested by CSS/test_submit_smoke_args.sh
# job_id and gen_njobs are how a run shards itself: task i generates the i-th of NJOBS
# slices. If that pair disagrees with the actual --array size, nothing errors -- the run
# exits 0 having quietly generated the wrong set of designs, which is the hardest kind of
# failure to notice because the output tree looks ordinary. Both directions are wrong:
#   * fewer array tasks than NJOBS  -> the missing shards are never generated;
#   * more array tasks than NJOBS   -> tasks collide on the same shard and redo each other's
#     work (with skip_existing on, the surplus tasks mostly idle through a GPU reservation);
#   * NJOBS > 1 with no --array at all -> only shard 0 runs, so (NJOBS-1)/NJOBS of the
#     designs are silently dropped. This is the one that bites, because the submission
#     looks complete and the job exits 0.
# Returns 1 and explains rather than calling exit, so the test can drive it directly.
check_array_njobs() {
    local n="${NJOBS:-1}"
    if [ -n "${SLURM_ARRAY_TASK_COUNT:-}" ] && [ "$SLURM_ARRAY_TASK_COUNT" != "$n" ]; then
        echo "FAIL: --array has $SLURM_ARRAY_TASK_COUNT task(s) but NJOBS=$n." >&2
        echo "      They must match: use --array=0-$((n - 1)), or set NJOBS=$SLURM_ARRAY_TASK_COUNT." >&2
        return 1
    fi
    if [ -z "${SLURM_ARRAY_TASK_ID:-}" ] && [ "$n" -gt 1 ]; then
        echo "FAIL: NJOBS=$n but this is not an array job, so only shard 0 would run and" >&2
        echo "      $((n - 1))/$n of the designs would never be generated -- with exit 0." >&2
        echo "      Submit with --array=0-$((n - 1)), or set NJOBS=1." >&2
        return 1
    fi
    return 0
}

# This script's own settings are ENVIRONMENT variables, read before sbatch ever runs the
# script. Everything after the script name is a Hydra override forwarded verbatim to
# rfd4-generate. Writing a setting in the second position looks identical to writing it in
# the first and is not:
#
#   CONFIG_NAME=x sbatch CSS/submit_generate.sh          # selects the arm
#   sbatch CSS/submit_generate.sh CONFIG_NAME=x          # Hydra: Could not override 'CONFIG_NAME'
#
# That is how the 10,000-design probe died (2026-11-19): the default arm was loaded, Hydra
# was handed a key no config defines, and the job failed after the environment was up. Hydra
# then suggests `+CONFIG_NAME=...`, which would APPEND the key and run the wrong arm quietly
# -- the worst available outcome, so this names the real fix instead.
check_env_in_args() {
    local arg key rc=0
    for arg in "$@"; do
        key="${arg%%=*}"
        case "$key" in
            CONFIG_NAME | NJOBS | SMOKE | REPO)
                echo "FAIL: '$arg' sets this script, not the model, so it must come BEFORE sbatch:" >&2
                echo "          $key=${arg#*=} sbatch [sbatch options] $0 [hydra overrides]" >&2
                echo "      Passed after the script name it reaches Hydra as a config key, which" >&2
                echo "      no config defines. Do NOT take Hydra's '+$key=...' suggestion: that" >&2
                echo "      appends the key and runs the DEFAULT arm without saying so." >&2
                rc=1
                ;;
        esac
    done
    return $rc
}
# <<< guards
check_array_njobs || exit 1
check_env_in_args "$@" || exit 1

cd "$REPO"
mkdir -p slurm_run_outputs
source .venv/bin/activate
export PROJECT_ROOT="$PWD"

# >>> timing report
# One record per array task, written on exit however the task ends. Two reasons it exists at
# all: Slurm's accounting database is the only other place this is kept, it is per-cluster, and
# `sacct -S` has to be aimed at the right window to see anything -- and the campaign's runs are
# split across Delta and Anvil, so there is no single query that answers "what did this arm
# cost". The record lands next to the job's .out file, which is already synced with the run.
#
# One file per task rather than appends to a shared CSV: ten to two hundred tasks appending
# concurrently on a parallel filesystem is exactly where interleaved short writes corrupt
# lines, and a lost timing line is worse than a directory with many small files in it.
#
# JSON Lines, so pooling is a cat and reading is one pandas call:
#   cat slurm_run_outputs/timings/*.json > timings.jsonl
#   python -c "import pandas;print(pandas.read_json('timings.jsonl',lines=True))"
TIMING_DIR="$PWD/slurm_run_outputs/timings"
mkdir -p "$TIMING_DIR"

# Bash substitution rather than a sed pipeline: the escaping needed to make sed escape
# backslashes is itself the kind of thing that survives review and then quietly emits nothing.
# $bs is one backslash, built with printf so this file carries no doubled escape of its own.
json_escape() {
    local s="${1:-}" bs
    bs="$(printf '%b' '\0134')"
    s="${s//$bs/$bs$bs}"
    s="${s//\"/$bs\"}"
    printf '%s' "$s"
}

timing_report() {
    local status="${1:-0}" t1 elapsed h m s hms file
    t1="$(date +%s)"
    elapsed=$((t1 - TIMING_T0))
    h=$((elapsed / 3600)); m=$(((elapsed % 3600) / 60)); s=$((elapsed % 60))
    hms="$(printf '%02d:%02d:%02d' "$h" "$m" "$s")"

    # %A_%a when this is an array task, the plain job id otherwise -- the same identity sacct
    # prints, so the two can be joined on it.
    local jobid="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}"
    [ -n "${SLURM_ARRAY_TASK_ID:-}" ] && jobid="${jobid}_${SLURM_ARRAY_TASK_ID}"

    echo
    echo "############ Wall clock ############"
    echo "elapsed     : $hms (${elapsed}s), exit status $status"
    echo "job         : $jobid on ${CLUSTER_LABEL} (${SLURM_JOB_PARTITION:-?})"

    file="$TIMING_DIR/${CONFIG_NAME}_${jobid}.json"
    {
        printf '{'
        printf '"cluster":"%s",'        "$(json_escape "$CLUSTER_LABEL")"
        printf '"config":"%s",'         "$(json_escape "$CONFIG_NAME")"
        printf '"job":"%s",'            "$(json_escape "$jobid")"
        printf '"array_task_id":"%s",'  "$(json_escape "${SLURM_ARRAY_TASK_ID:-}")"
        printf '"njobs":%s,'            "${NJOBS:-1}"
        printf '"partition":"%s",'      "$(json_escape "${SLURM_JOB_PARTITION:-}")"
        printf '"account":"%s",'        "$(json_escape "${SLURM_JOB_ACCOUNT:-}")"
        printf '"node":"%s",'           "$(json_escape "${SLURM_JOB_NODELIST:-}")"
        printf '"gpus":"%s",'           "$(json_escape "${SLURM_GPUS_ON_NODE:-${SLURM_JOB_GPUS:-}}")"
        printf '"smoke":%s,'            "$(case "${SMOKE:-0}" in 1 | true | TRUE | yes | YES) echo true ;; *) echo false ;; esac)"
        printf '"out_dir":"%s",'        "$(json_escape "${SMOKE_OUT_DIR:-${CONFIG_OUT_DIR:-}}")"
        printf '"overrides":"%s",'      "$(json_escape "$TIMING_ARGS")"
        printf '"start":"%s",'          "$(json_escape "$TIMING_START_ISO")"
        printf '"end":"%s",'            "$(json_escape "$(date -Is)")"
        printf '"elapsed_s":%s,'        "$elapsed"
        printf '"elapsed":"%s",'        "$hms"
        printf '"gpu_hours":%s.%03d,'   "$((elapsed / 3600))" "$(((elapsed % 3600) * 1000 / 3600))"
        printf '"exit_status":%s'       "$status"
        printf '}\n'
    } > "$file"
    echo "recorded    : $file"
    :
}

# `$?` is expanded when the trap fires, not now, so the record carries the real exit status --
# including the guard failures and the smoke gate's exit 1, which are the runs whose cost is
# easiest to lose track of.
TIMING_ARGS="$*"
trap 'timing_report $?' EXIT
# <<< timing report

# generate.py calls load_dotenv(), which reads .env.<CLUSTER> and falls back to .env.local.
# Leaving CLUSTER unset selects .env.local -- the same profile interactive runs use. Setting
# it to a cluster name would demand .env.<name>, which is how the first batch attempt failed.
unset CLUSTER || true

# The whole point: with this unset, the oracle calls apptainer directly instead of dialling a
# daemon socket that does not exist here.
unset AF3_LAUNCHER_SOCKET_TEMPLATE || true

# Fail here rather than after generation finishes. The AF3 paths live in .env.local; sourcing
# it now only to check them is harmless, since python-dotenv does not overwrite values already
# in the environment and these are the same ones it would set.
set -a; source .env.local; set +a
for var in AF3_IMAGE AF3_MODEL_DIR AF3_JAX_CACHE_DIR; do
    if [ -z "${!var:-}" ]; then
        echo "FAIL: $var is empty. Refolding would be skipped and the run would produce" >&2
        echo "      designs with no oracle metrics. Set it in .env.local." >&2
        exit 1
    fi
done
if [ ! -e "$AF3_IMAGE" ]; then
    echo "FAIL: AF3_IMAGE does not exist: $AF3_IMAGE" >&2
    exit 1
fi

# Fail on a bad arm name here, not after Hydra has loaded a checkpoint.
CONFIG_FILE="$PWD/CSS/configs/${CONFIG_NAME}.yaml"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "FAIL: no such config: $CONFIG_FILE" >&2
    echo "      Available:" >&2
    ls "$PWD/CSS/configs" >&2
    exit 1
fi

# Best effort, for the timing record only: the arm's own out_dir, so a recorded run can be
# matched to the tree it wrote. Empty when the key comes from a composed default rather than
# this file, and a SMOKE run's fresh directory wins over it.
CONFIG_OUT_DIR="$(sed -n 's/^out_dir:[[:space:]]*//p' "$CONFIG_FILE" | head -1 \
    | sed "s|\${oc.env:PROJECT_ROOT}|$PROJECT_ROOT|g; s/^\"//; s/\"$//; s/^'//; s/'$//")"
# The hotspot-sampling arms instantiate a transform per replicate. Two ways that costs a
# whole allocation if it is wrong, both caught in seconds here:
#   * KeepSeedHotspots lives under src/rfproteina/transforms/sampling/, which .gitignore:9
#     excludes. The module is force-added so a `git pull` brings it, but the directory is
#     still ignored, and a checkout without it kills the add arm on import at the first batch;
#   * SampleNearbyHotspots raises rather than reducing the count when the patch cannot fit
#     the requested hotspots on these inputs.
# This also reports whether the seeds are kept or replaced, which is the difference between
# the add and swap arms and is otherwise invisible in the run's own output.
if grep -q "sampling_transforms" "$CONFIG_FILE"; then
    echo "############ sampling arm: hotspot preflight ############"
    python CSS/preflight_hotspot_patch.py 10 --config "$CONFIG_NAME"
    echo
fi

# >>> smoke args -- extracted and unit-tested by CSS/test_submit_smoke_args.sh
# SMOKE=1 turns this into stage 1 of CSS/run_smoke.sh as a batch job: 3 replicates over the
# 4 inputs, batch size 1, a fresh out_dir, and the column audit afterwards. Without it the
# arms run at their configured size -- 100 replicates x 4 inputs = 400 designs and 800 AF3
# refolds, which is a production run, not a check that the plumbing works.
#
# Anything you pass explicitly wins: each default is added only if that key is absent from
# the arguments, because Hydra is not relied on to resolve a duplicated override.
smoke_enabled() { case "${SMOKE:-0}" in 1 | true | TRUE | yes | YES) return 0 ;; *) return 1 ;; esac; }

has_override() {
    local key="$1" arg
    shift
    for arg in "$@"; do
        [ "${arg%%=*}" = "$key" ] && return 0
    done
    return 1
}

# How many structure files does this config's `inputs` point at? The smoke expectation is
# per input, and the arms do not agree on that number: the four-framework configs name a
# DIRECTORY of 4 CIFs, while the 1000-design depth probes name a SINGLE vobarilizumab CIF.
# A hardcoded 12 fails a perfectly good single-input smoke run, which is the same false
# alarm the gate exists to prevent in the other direction. Echoes 0 when the config cannot
# be read or the path does not resolve, and the caller then only requires "something".
count_inputs() {
    local file="${CONFIG_FILE:-}" spec root n=0 p
    if [ -z "$file" ] || [ ! -f "$file" ]; then echo 0; return; fi
    root="${PROJECT_ROOT:-$PWD}"
    spec=$(sed -n 's/^inputs:[[:space:]]*//p' "$file" | head -1 | sed "s|\${oc.env:PROJECT_ROOT}|$root|g")
    spec=$(printf '%s' "$spec" | sed 's/^"//; s/"$//; s/^'\''//; s/'\''$//')
    if [ -z "$spec" ]; then echo 0; return; fi
    if [ -d "$spec" ]; then
        find "$spec" -maxdepth 1 -type f \( -name '*.cif' -o -name '*.cif.zst' \) | wc -l
    elif [ -f "$spec" ]; then
        echo 1
    else
        # A glob, or a path this cannot resolve from here.
        for p in $spec; do [ -e "$p" ] && n=$((n + 1)); done
        echo "$n"
    fi
}

build_smoke_args() {
    # Sets SMOKE_ARGS, SMOKE_OUT_DIR and SMOKE_EXPECTED from the caller's arguments.
    SMOKE_ARGS=()
    SMOKE_OUT_DIR=""
    SMOKE_EXPECTED=0
    local arg
    if has_override num_replicates "$@"; then
        # Their number, so this script cannot say how many designs to expect.
        :
    else
        SMOKE_ARGS+=(num_replicates=3)
        # 3 replicates per input: 12 for the four-framework arms, 3 for the single-input
        # depth probes. Each design writes a design plus two refolds, so counting one file
        # per design keeps this a lower bound rather than an exact figure.
        local n_inputs
        n_inputs=$(count_inputs)
        if [ "$n_inputs" -gt 0 ]; then
            SMOKE_EXPECTED=$((3 * n_inputs))
        else
            SMOKE_EXPECTED=1
        fi
        # Under --array the designs are split across tasks, so one task's share is smaller
        # and the full expectation would fail every task. Only "something" holds then.
        [ "${NJOBS:-1}" -gt 1 ] && SMOKE_EXPECTED=1
    fi
    has_override batch_size "$@" || SMOKE_ARGS+=(batch_size=1)
    has_override seed "$@" || SMOKE_ARGS+=(seed=42)
    has_override verbose "$@" || SMOKE_ARGS+=(verbose=true)
    has_override prevalidate_inputs "$@" || SMOKE_ARGS+=(prevalidate_inputs=true)
    if has_override out_dir "$@"; then
        for arg in "$@"; do
            [ "${arg%%=*}" = out_dir ] && SMOKE_OUT_DIR="${arg#*=}"
        done
    else
        # Fresh every time. `skip_existing` defaults to true, so a reused directory keeps the
        # previous run's structures and the audit below would read stale results.
        SMOKE_OUT_DIR="$PWD/runs/smoke_batch_${CONFIG_NAME}_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
        SMOKE_ARGS+=("out_dir=$SMOKE_OUT_DIR")
    fi
}
# <<< smoke args

SMOKE_ARGS=()
SMOKE_OUT_DIR=""
SMOKE_EXPECTED=0
if smoke_enabled; then
    build_smoke_args "$@"
fi

echo "repo        : $PWD"
echo "config      : $CONFIG_NAME"
echo "array task  : ${SLURM_ARRAY_TASK_ID:-0} of $NJOBS"
echo "AF3 image   : $AF3_IMAGE"
echo "AF3 daemon  : disabled (direct apptainer path)"
if smoke_enabled; then
    echo "size        : SMOKE -- ${SMOKE_ARGS[*]}"
else
    echo "size        : the config's own (production; SMOKE=1 for a 12-design check)"
fi
echo "overrides   : $*"

# `${a[@]+"${a[@]}"}` because an empty array under `set -u` is an error on bash 4.2, which is
# what the compute nodes have.
rfd4-generate \
    --config-path "$PWD/CSS/configs" \
    --config-name "$CONFIG_NAME" \
    save_chain_of_interest=true \
    job_id="${SLURM_ARRAY_TASK_ID:-0}" \
    gen_njobs="$NJOBS" \
    ${SMOKE_ARGS[@]+"${SMOKE_ARGS[@]}"} \
    "$@"

if smoke_enabled && [ -n "$SMOKE_OUT_DIR" ]; then
    echo
    echo "############ Smoke gate: did it really produce structures? ############"
    # Exiting 0 having produced nothing is this campaign's signature failure -- a bad input
    # path or an input glob matching no file does exactly that. The audit below would not
    # notice, since it reads whatever rows exist.
    n=$(find "$SMOKE_OUT_DIR" -type f \( -name '*.cif' -o -name '*.cif.zst' \) 2>/dev/null | wc -l)
    if [ "$SMOKE_EXPECTED" -gt 0 ] && [ "$n" -lt "$SMOKE_EXPECTED" ]; then
        echo "FAIL: produced $n structure file(s), expected at least $SMOKE_EXPECTED" >&2
        echo "      under $SMOKE_OUT_DIR. The job exited 0, so this would have looked fine." >&2
        exit 1
    fi
    if [ "$n" -lt 1 ]; then
        echo "FAIL: produced no structure files under $SMOKE_OUT_DIR" >&2
        exit 1
    fi
    echo "ok: $n structure file(s) under $SMOKE_OUT_DIR"

    echo
    echo "############ Audit: are the metric columns there? ############"
    python CSS/audit_smoke_metrics.py "$SMOKE_OUT_DIR"
fi
