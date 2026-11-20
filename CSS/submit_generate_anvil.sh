#!/bin/bash
# Submit a CSS generation run on Anvil (Purdue / ACCESS).
#
# Why this is a wrapper and not a copy
# ------------------------------------
# Everything about HOW a run executes -- the AF3 daemon bypass, the AF3_* guards, the
# sampling-arm hotspot preflight, the SMOKE sizing with its post-run structure-count gate and
# metric audit -- lives in CSS/submit_generate.sh and is cluster-neutral. The only
# Delta-specific parts of that file are its #SBATCH header and the REPO default.
#
# SLURM reads directives only from the file handed to sbatch, so this file supplies Anvil's
# header and then hands off. The #SBATCH lines inside submit_generate.sh are ordinary comments
# when bash runs it as a sub-script, so they are inert here. Keeping this a wrapper means a fix
# to the smoke logic lands on both clusters at once instead of drifting apart -- the same reason
# the two clusters share one branch rather than having an `anvil` branch of their own.
#
# Usage, from the repository root. Extra arguments pass through to rfd4-generate as Hydra
# overrides, exactly as for the Delta script:
#
#   mkdir -p slurm_run_outputs                               # once, before the first submit
#   sbatch --array=0-9 CSS/submit_generate_anvil.sh          # the groove arm, as configured
#   SMOKE=1 sbatch --time=01:00:00 CSS/submit_generate_anvil.sh
#   CONFIG_NAME=css_vobarilizumab_1000_generate sbatch --array=0-9 CSS/submit_generate_anvil.sh
#
# Default arm: the four-hotspot groove-spanning epitope
# E165.OE2 / L204.CD2 / L206.CD2 / Q241.NE2, from
# CSS/configs/css_groove_E165OE2_L204CD2_L206CD2_Q241NE2_generate.yaml -- 250 replicates across
# the 4 VHH frameworks in CSS/inputs/cifs_E165OE2_L204CD2_L206CD2_Q241NE2/ = 1000 designs and
# 2000 AF3 refolds, medium / front-loaded recycling, batch_size 2, compressed outputs.
#
# Sizing: the config sets gen_njobs 10, so a production submission needs --array=0-9 and the
# NJOBS default of 10. The guards below refuse a mismatch rather than let it through: without
# an array, job_id stays 0 and nine tenths of the designs are silently never generated, which
# exits 0 and looks like a completed run. SMOKE=1 defaults NJOBS to 1 instead, since a smoke
# run is 3 replicates per input (12 designs here) in a single task.
#
# Account and partition: bio260655-ai maps to the `ai` partition (h[000-020]). If you also hold
# bio260655-gpu, the `gpu` partition's A100s are the same GPU generation as Delta, which is the
# better choice for anything you intend to compare against runs/css_high:
#
#   sbatch -A bio260655-gpu -p gpu --array=0-9 CSS/submit_generate_anvil.sh
#
# Resuming, exactly as on Delta: resubmit with the SAME out_dir. `skip_existing` defaults to
# true, so finished designs are skipped and refolds with an existing confidence JSON are not
# repeated. On a congested queue a short --time plus a resubmit backfills far better than one
# long reservation.

#SBATCH --job-name=css-groove
#SBATCH --account=bio260655-ai
#SBATCH --partition=ai
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --time=04:00:00
#SBATCH --signal=SIGUSR1@600
#SBATCH --output=slurm_run_outputs/css_groove_%A_%a.out
#SBATCH --error=slurm_run_outputs/css_groove_%A_%a.err

set -euo pipefail

export REPO="${REPO:-/home/x-agurajapu/RFD4-Proteina}"
export CONFIG_NAME="${CONFIG_NAME:-css_groove_E165OE2_L204CD2_L206CD2_Q241NE2_generate}"

smoke_enabled() { case "${SMOKE:-0}" in 1 | true | TRUE | yes | YES) return 0 ;; *) return 1 ;; esac; }
if smoke_enabled; then
    export NJOBS="${NJOBS:-1}"
else
    export NJOBS="${NJOBS:-10}"
fi

# Catch an array/gen_njobs mismatch here rather than discovering it in the design count. Both
# directions are silent failures: too few tasks drops designs, too many re-runs the same shard.
if [ -n "${SLURM_ARRAY_TASK_COUNT:-}" ] && [ "$SLURM_ARRAY_TASK_COUNT" != "$NJOBS" ]; then
    echo "FAIL: --array has $SLURM_ARRAY_TASK_COUNT task(s) but NJOBS=$NJOBS." >&2
    echo "      They must match: use --array=0-$((NJOBS - 1)), or set NJOBS=$SLURM_ARRAY_TASK_COUNT." >&2
    exit 1
fi
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ] && [ "$NJOBS" -gt 1 ]; then
    echo "FAIL: NJOBS=$NJOBS but this is not an array job, so only shard 0 would run and" >&2
    echo "      $((NJOBS - 1))/$NJOBS of the designs would never be generated -- with exit 0." >&2
    echo "      Submit with --array=0-$((NJOBS - 1)), or set NJOBS=1." >&2
    exit 1
fi

echo "cluster     : anvil (${SLURM_JOB_PARTITION:-?}, account ${SLURM_JOB_ACCOUNT:-?})"
exec bash "$REPO/CSS/submit_generate.sh" "$@"
