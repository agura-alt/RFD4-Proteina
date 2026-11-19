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

REPO="${REPO:-/u/agurajapu/RFD4-Proteina}"
CONFIG_NAME="${CONFIG_NAME:-css_cdr123_L204CD2_L206CD2_Y246CE1_high_generate}"
NJOBS="${NJOBS:-1}"   # must equal the --array size when fanning out

cd "$REPO"
mkdir -p slurm_run_outputs
source .venv/bin/activate
export PROJECT_ROOT="$PWD"

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

echo "repo        : $PWD"
echo "config      : $CONFIG_NAME"
echo "array task  : ${SLURM_ARRAY_TASK_ID:-0} of $NJOBS"
echo "AF3 image   : $AF3_IMAGE"
echo "AF3 daemon  : disabled (direct apptainer path)"
echo "overrides   : $*"

rfd4-generate \
    --config-path "$PWD/CSS/configs" \
    --config-name "$CONFIG_NAME" \
    save_chain_of_interest=true \
    job_id="${SLURM_ARRAY_TASK_ID:-0}" \
    gen_njobs="$NJOBS" \
    "$@"
