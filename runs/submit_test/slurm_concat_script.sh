#!/bin/bash
#SBATCH --job-name concat-inference
#SBATCH --account bjua-delta-gpu
#SBATCH --partition gpuA100x4
#SBATCH --time 04:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 16g
#SBATCH --output ./slurm_run_outputs/slurm_concat.out
#SBATCH --error ./slurm_run_outputs/slurm_concat.err
#SBATCH --open-mode append

export PYTHONPATH=/u/agurajapu/RFD4-Proteina:$PYTHONPATH
export PYTHONFAULTHANDLER=1

echo "Running result concatenation for: inference"
echo "Concat dir: /u/agurajapu/RFD4-Proteina/runs/submit_test/inference_metrics"
# Cluster identity — read by exec_core.sh (shebang) to source .ipd/.env.${CLUSTER}
export CLUSTER=delta
export CONTAINER_TYPE=none

"/u/agurajapu/RFD4-Proteina/src/rfproteina/inference_utils/aggregate_metrics.py" --input "/u/agurajapu/RFD4-Proteina/runs/submit_test/inference_metrics" --cluster-with-foldseek
echo "Concatenation complete!"
