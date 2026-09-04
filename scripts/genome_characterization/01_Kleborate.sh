#!/bin/bash
#SBATCH --job-name=GPSKp_Kleborate
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
CONTAINER="${KP_KLEBORATE_CONTAINER:-$PWD/containers/kleborate.sif}"

RESULTS="${BASE}/kleborate/results"
STATUS="${BASE}/kleborate/status"

mkdir -p "$RESULTS" "$STATUS"

TASK_ID="${SLURM_ARRAY_TASK_ID}"
LINE=$((TASK_ID + 1))

SRR=$(sed -n "${LINE}p" "$MANIFEST" | cut -f1)
ASSEMBLY=$(sed -n "${LINE}p" "$MANIFEST" | cut -f2)

if [[ -z "$SRR" || -z "$ASSEMBLY" ]]; then
    echo "ERROR: Could not read manifest entry for task ${TASK_ID}"
    exit 1
fi

OUTDIR="${RESULTS}/${SRR}"

echo "=================================================="
echo "GPS-Kp Kleborate"
echo "Sample: ${SRR}"
echo "Assembly: ${ASSEMBLY}"
echo "Output directory: ${OUTDIR}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "=================================================="

if [[ ! -s "$ASSEMBLY" ]]; then
    echo "ERROR: Assembly is missing or empty:"
    echo "$ASSEMBLY"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

if [[ ! -s "$CONTAINER" ]]; then
    echo "ERROR: Kleborate container is missing or empty:"
    echo "$CONTAINER"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

if [[ -d "$OUTDIR" ]] && find "$OUTDIR" -type f -size +0c | grep -q .; then
    echo "Existing non-empty result directory found; skipping."
    rm -f "${STATUS}/${SRR}.failed"
    touch "${STATUS}/${SRR}.complete"
    exit 0
fi

rm -rf "$OUTDIR"

apptainer exec \
    "$CONTAINER" \
    kleborate \
    --assemblies "$ASSEMBLY" \
    --outdir "$OUTDIR" \
    --preset kpsc

if [[ ! -d "$OUTDIR" ]] || ! find "$OUTDIR" -type f -size +0c | grep -q .; then
    echo "ERROR: Kleborate did not create valid results for ${SRR}"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

rm -f "${STATUS}/${SRR}.failed"
touch "${STATUS}/${SRR}.complete"

echo "=================================================="
echo "Completed Kleborate successfully for ${SRR}"
echo "Finished: $(date)"
echo "=================================================="


