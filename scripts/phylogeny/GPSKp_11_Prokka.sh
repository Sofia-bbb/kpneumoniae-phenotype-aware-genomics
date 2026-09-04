#!/bin/bash
#SBATCH --job-name=GPSKp_Prokka
#SBATCH --account=def-walshj
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/Prokka_%A_%a.out
#SBATCH --error=logs/Prokka_%A_%a.err

set -euo pipefail

BASE="${KP_WORKDIR:-$PWD/work}"
MANIFEST="${KP_PHYLOGENY_MANIFEST:-$PWD/config/assemblies_phylogeny.tsv}"
OUTROOT="${KP_PROKKA_OUTROOT:-${BASE}/prokka/annotations}"
STATUS="${KP_PROKKA_STATUS:-${BASE}/prokka/status}"

mkdir -p "${OUTROOT}" "${STATUS}"

module --force purge
module load StdEnv/2020
module load gcc/9.3.0
module load prokka/1.14.5

# Manifest has a header; array task 1 reads line 2.
LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 1))

LINE=$(sed -n "${LINE_NUMBER}p" "${MANIFEST}")

if [[ -z "${LINE}" ]]; then
    echo "ERROR: No manifest entry for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

SRR=$(printf '%s\n' "${LINE}" | cut -f1)
ASSEMBLY=$(printf '%s\n' "${LINE}" | cut -f2)

if [[ -z "${SRR}" || -z "${ASSEMBLY}" ]]; then
    echo "ERROR: Could not parse SRR or assembly path from line ${LINE_NUMBER}" >&2
    exit 1
fi

if [[ ! -s "${ASSEMBLY}" ]]; then
    echo "ERROR: Assembly missing or empty: ${ASSEMBLY}" >&2
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

OUTDIR="${OUTROOT}/${SRR}"

rm -f \
    "${STATUS}/${SRR}.complete" \
    "${STATUS}/${SRR}.failed"

# Skip a sample that already has a valid GFF output.
if [[ -s "${OUTDIR}/${SRR}.gff" ]]; then
    echo "Existing valid annotation found for ${SRR}; skipping."
    touch "${STATUS}/${SRR}.complete"
    exit 0
fi

rm -rf "${OUTDIR}"

echo "============================================================"
echo "Sample: ${SRR}"
echo "Assembly: ${ASSEMBLY}"
echo "Output: ${OUTDIR}"
echo "Started: $(date)"
echo "============================================================"

if prokka \
    --outdir "${OUTDIR}" \
    --prefix "${SRR}" \
    --locustag "${SRR}" \
    --kingdom Bacteria \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --force \
    --quiet \
    "${ASSEMBLY}"
then
    if [[ -s "${OUTDIR}/${SRR}.gff" ]]; then
        touch "${STATUS}/${SRR}.complete"
        echo "Completed successfully: ${SRR}"
    else
        touch "${STATUS}/${SRR}.failed"
        echo "ERROR: Prokka finished but GFF is missing or empty." >&2
        exit 1
    fi
else
    touch "${STATUS}/${SRR}.failed"
    echo "ERROR: Prokka failed for ${SRR}" >&2
    exit 1
fi

echo "Finished: $(date)"


