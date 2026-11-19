#!/bin/bash
# Build and score the native 1G8 complex as a positive control for this campaign's gates.
#
# Why
# ---
# Two runs, 112 designs, two recycling budgets, and oracle_ppi.min_ipae_coi has not moved:
# 15.8-24.5 at 180 trunk passes, 20.1 +/- 3.4 at 300, against a gate of 1.2. At that spread the
# gate sits ~5.6 sigma below the mean, so scaling does not reach it -- the expected best of
# 10,000 designs is about 7, against an observed best of 8.96 in 100. The distribution would
# have to move, not be sampled deeper.
#
# That leaves one question unanswered: is 1.2 reachable on this antigen at all? 1G8 is a
# nanobody with a solved structure against this target, and it contacts all three of this
# campaign's conditioned hotspot atoms (3.02 / 3.55 / 3.61 A, inside the gate's own 4.5 A). It
# is a binder that demonstrably works, at this epitope. Scoring it through the same oracle
# gives the ceiling.
#
#   native near 1    the gate is reachable; the designs genuinely do not dock.
#   native near 20   a crystallographic binder fails the gate, so min_ipae cannot separate
#                    binders from non-binders here and every 0/112 in this campaign is a
#                    statement about the metric. Switch to a filter that varies with design
#                    quality -- apo/holo has a median of 2.04 A and a quartile under 0.81 A,
#                    where min_ipae is pinned near 20 regardless.
#
# Two AF3 calls. Submitted as a batch job rather than run interactively because an srun --pty
# waiting in the queue dies with its SSH session, which has already cost one attempt.
#
# Usage, from the repository root. Extra arguments pass through to rfd4-evaluate as Hydra
# overrides:
#   sbatch CSS/submit_native_control.sh
#   sbatch --time=02:00:00 CSS/submit_native_control.sh

#SBATCH --job-name=native-control
#SBATCH --account=bjua-delta-gpu
#SBATCH --partition=gpuA100x4
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --time=01:00:00
#SBATCH --output=slurm_run_outputs/native_control_%j.out
#SBATCH --error=slurm_run_outputs/native_control_%j.err

set -euo pipefail

REPO="${REPO:-/u/agurajapu/RFD4-Proteina}"
CONFIG_NAME="${CONFIG_NAME:-css_native_control_evaluate}"
OUT_DIR="${OUT_DIR:-$REPO/runs/native_control_annotated}"

cd "$REPO"
mkdir -p slurm_run_outputs
source .venv/bin/activate
export PROJECT_ROOT="$PWD"

# Leaving CLUSTER unset makes load_dotenv() read .env.local, the profile interactive runs use.
unset CLUSTER || true
# rfd4-evaluate drives the same AF3Oracle as generation. With this unset the oracle calls
# apptainer directly instead of dialling a daemon socket, which is the path that works here --
# the launcher's pre-spawned daemon never came up on this cluster (job 24147097).
unset AF3_LAUNCHER_SOCKET_TEMPLATE || true

# Fail before spending the allocation. An unset AF3 environment does not raise on its own:
# the run would finish having scored nothing, and report success.
set -a; source .env.local; set +a
for var in AF3_IMAGE AF3_MODEL_DIR AF3_JAX_CACHE_DIR; do
    if [ -z "${!var:-}" ]; then
        echo "FAIL: $var is empty; set it in .env.local" >&2
        exit 1
    fi
done
[ -e "$AF3_IMAGE" ] || { echo "FAIL: AF3_IMAGE does not exist: $AF3_IMAGE" >&2; exit 1; }

echo "=============================== build ==============================="
# Prints its CDR selection, the residue count and the hotspot contact report, and raises if
# renumbering would merge 1G8's antibody-numbering insertions (52A, 82A-C, 100A-D) -- which
# silently dropped the chain from 124 residues to 116 and made AF3 fold it to nonsense
# (pLDDT 39.6) on the first attempt. Expect "renumbered 124 residues 1..124".
python CSS/build_native_control_annotated.py

echo "=============================== score ==============================="
rfd4-evaluate --config-path "$PWD/CSS/configs" --config-name "$CONFIG_NAME" "$@"

echo "============================== results =============================="
python - "$OUT_DIR" <<'PY'
import sys, glob
import pandas as pd

csvs = sorted(glob.glob(f"{sys.argv[1]}/evaluation_metrics/per_example_metrics*.csv"))
if not csvs:
    print(f"no per_example_metrics*.csv under {sys.argv[1]}/evaluation_metrics")
    raise SystemExit(1)
d = pd.concat([pd.read_csv(c) for c in csvs], ignore_index=True)

# oracle_apo first: if AF3 could not fold the native nanobody alone, nothing else is
# interpretable. A canonical VHH should be well above 80 here; 31 was the mangled-input run.
order = [
    "oracle_apo.apo.ptm_coi",
    "oracle_apo.apo.plddt_coi",
    "oracle_ppi.ptm_coi",
    "oracle_ppi.plddt_coi",
    "oracle_ppi.iptm_coi",
    "oracle_ppi.min_ipae_coi",
    "refolded_apo_holo_rmsd.apo.rmsd",
    "refolded_apo_holo_cdr_rmsd.apo.rmsd",
    "hotspot_obedience.hotspot_obedience",
    "min_ipae.all_metrics_present",
    "min_ipae.pass",
    "hotspot_adherence.all_metrics_present",
    "hotspot_adherence.pass",
    "apo_holo_consistency.all_metrics_present",
    "apo_holo_consistency.pass",
]
present = [c for c in order if c in d.columns]
print(d[present].T.to_string())
missing = [c for c in order if c not in d.columns]
if missing:
    print("\nMISSING COLUMNS:", ", ".join(missing))

print("\n--- read in this order ---")
print("1. oracle_apo.apo.ptm_coi  -- below ~70 means AF3 did not fold the native binder;")
print("                              stop, nothing downstream means anything.")
print("2. *.all_metrics_present   -- any 0 makes that criterion indeterminate, not failed.")
print("3. oracle_ppi.min_ipae_coi -- vs the campaign's 20.1 +/- 3.4 and the 1.2 gate.")
print("4. refolded_apo_holo_*     -- vs the campaign median 2.04 A / 2.0 A criterion.")
PY
