#!/bin/bash
#SBATCH --job-name inference
#SBATCH --time 04:00:00
#SBATCH --nodes 1
#SBATCH --gpus-per-node 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 8
#SBATCH --open-mode append
#SBATCH --requeue
#SBATCH --output ./slurm_run_outputs/slurm.out
#SBATCH --error ./slurm_run_outputs/slurm.err

# Raise the open-file soft limit to the hard cap on hosts with low defaults.
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
cd /u/agurajapu/RFD4-Proteina

# Environment setup (all clusters)
export PYTHONPATH=/u/agurajapu/RFD4-Proteina:$PYTHONPATH
export PYTHONFAULTHANDLER=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

# Job array configuration
ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
GPUS_PER_JOB=1
NJOBS=1

# --exclusive below gives each step its own GPU, but caps it at one CPU unless asked.
CPUS_PER_STEP=$(( ${SLURM_CPUS_ON_NODE:-8} / GPUS_PER_JOB ))
CPUS_PER_STEP=$(( CPUS_PER_STEP >= 2 ? CPUS_PER_STEP : 2 ))  # rule of thumb: never < 2 CPUs/GPU

echo "Array task: $ARRAY_TASK_ID, GPUs per job: $GPUS_PER_JOB, total jobs: $NJOBS, CPUs per step: $CPUS_PER_STEP"

# Cluster identity — read by exec_core.sh (shebang) to source .ipd/.env.${CLUSTER}
export CLUSTER=local
export CONTAINER_TYPE=none

export MASTER_PORT=$((1024 + RANDOM % 64512))
export MASTER_ADDR=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n 1)

# Verify scontrol works (needed for auto-requeue)
scontrol show job $SLURM_JOB_ID > /dev/null 2>&1 && echo '[slurm] scontrol works - auto-requeue enabled' || echo '[slurm] WARNING: scontrol not working - auto-requeue disabled'

# Pre-launch AF3 daemons (one supervisor per node, sockets in each node's /tmp).
# Foreground srun + shell background + --overlap keeps the supervisor alive without slurm killing it.
export AF3_IMAGE="/u/agurajapu/RFD4-Proteina/data_path/openfold3/runtime/of3_af3_refolding_env-cuda13.sif"
export AF3_JAX_CACHE_DIR="/u/agurajapu/RFD4-Proteina/data_path/openfold3/jax_cache"
export AF3_MODEL_DIR="/u/agurajapu/RFD4-Proteina/data_path/af3_models"
AF3_OUT=${RUN_DIR:-.}/slurm_run_outputs && mkdir -p "$AF3_OUT"
srun --ntasks-per-node=1 --nodes=$SLURM_JOB_NUM_NODES --overlap --kill-on-bad-exit=0 \
  /usr/bin/python3 /u/agurajapu/RFD4-Proteina/src/rfproteina/af3/af3_daemon.py launch \
  --num-gpus ${SLURM_GPUS_PER_NODE:-1} --log-dir "$AF3_OUT" \
  --skip-if-no-env 2>>"$AF3_OUT/af3_daemons.log" >/dev/null &
AF3_DAEMON_SRUN_PID=$!
trap 'kill -TERM $AF3_DAEMON_SRUN_PID 2>/dev/null; wait $AF3_DAEMON_SRUN_PID 2>/dev/null' EXIT
export AF3_LAUNCHER_SOCKET_TEMPLATE="/tmp/af3_${SLURM_JOB_ID}/af3_daemon_rank{rank}.sock"

# Run generation - each GPU gets its own job_id slice of the dataset
GEN_PIDS=()
for LOCAL_ID in $(seq 0 $((GPUS_PER_JOB - 1))); do
    JOB_ID=$((ARRAY_TASK_ID * GPUS_PER_JOB + LOCAL_ID))
    if [ $JOB_ID -lt $NJOBS ]; then
        echo "Running generation with job_id=$JOB_ID on GPU $LOCAL_ID"
        AF3_DAEMON_RANK=$LOCAL_ID AF2_DAEMON_RANK=$LOCAL_ID srun --exclusive --ntasks=1 --gpus=1 \
            --cpus-per-task=$CPUS_PER_STEP --cpu-bind=none \
            --kill-on-bad-exit /u/agurajapu/RFD4-Proteina/src/rfproteina/generate.py \
            --config-path /u/agurajapu/RFD4-Proteina/CSS/configs --config-name css_cdr123_L204CD2_L206CD2_Y246CE1_high_generate \
             save_chain_of_interest=true \
            out_dir=/u/agurajapu/RFD4-Proteina/runs/css_cdr123_L204CD2_L206CD2_Y246CE1_high \
            hardware.ngpus_per_node_=1 \
            hardware.nnodes_=1 \
            job_id=$JOB_ID \
            gen_njobs=$NJOBS &
        GEN_PIDS+=($!)
    fi
done
wait "${GEN_PIDS[@]}"
