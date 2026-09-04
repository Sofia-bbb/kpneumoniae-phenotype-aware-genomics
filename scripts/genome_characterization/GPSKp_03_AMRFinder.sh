#!/bin/bash
#SBATCH --job-name=GPSKp_AMRFinder
#SBATCH --account=def-walshj
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

module purge
module load StdEnv/2023
module load apptainer/1.4.5

BASE="${KP_WORKDIR:-$PWD/work}"
MANIFEST="${KP_ASSEMBLY_MANIFEST:-$PWD/config/assemblies.tsv}"
CONTAINER="${KP_AMRFINDER_CONTAINER:-$PWD/containers/kleborate.sif}"
RESULTS="${BASE}/amrfinder/results"
STATUS="${BASE}/amrfinder/status"

mkdir -p "$RESULTS" "$STATUS"

TASK_ID="${SLURM_ARRAY_TASK_ID}"
LINE=$((TASK_ID + 1))

SRR=$(sed -n "${LINE}p" "$MANIFEST" | cut -f1)
ASSEMBLY=$(sed -n "${LINE}p" "$MANIFEST" | cut -f2)

if [[ -z "$SRR" || -z "$ASSEMBLY" ]]; then
    echo "ERROR: Could not read manifest entry for task ${TASK_ID}"
    exit 1
fi

OUT="${RESULTS}/${SRR}_amrfinder.tsv"

echo "=================================================="
echo "GPS-Kp AMRFinderPlus"
echo "Sample: ${SRR}"
echo "Assembly: ${ASSEMBLY}"
echo "Output: ${OUT}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "=================================================="

if [[ ! -s "$ASSEMBLY" ]]; then
    echo "ERROR: Assembly is missing or empty: $ASSEMBLY"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

if [[ ! -s "$CONTAINER" ]]; then
    echo "ERROR: Container is missing or empty: $CONTAINER"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

if [[ -s "$OUT" ]]; then
    echo "Existing non-empty AMRFinderPlus result found; skipping."
    rm -f "${STATUS}/${SRR}.failed"
    touch "${STATUS}/${SRR}.complete"
    exit 0
fi

rm -f "$OUT"

apptainer exec \
    "$CONTAINER" \
    amrfinder \
    -n "$ASSEMBLY" \
    -O Klebsiella_pneumoniae \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -o "$OUT"

if [[ ! -s "$OUT" ]]; then
    echo "ERROR: AMRFinderPlus result was not created for ${SRR}"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

rm -f "${STATUS}/${SRR}.failed"
touch "${STATUS}/${SRR}.complete"

echo "=================================================="
echo "Completed AMRFinderPlus successfully for ${SRR}"
echo "Finished: $(date)"
echo "=================================================="
