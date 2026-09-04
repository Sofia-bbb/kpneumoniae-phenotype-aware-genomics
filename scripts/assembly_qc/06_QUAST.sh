#!/bin/bash
#SBATCH --job-name=KP_QUAST
#SBATCH --account=def-walshj
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-5%2
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# ============================================================
# Script: 06_QUAST.sh
# Project: Klebsiella pneumoniae genotype-phenotype study
# Authors: Sofia B. Mohamed and Jacob P. Walsh
#
# Purpose:
#   Assess the quality of individual SPAdes assemblies using
#   QUAST.
#
# Input:
#   work/assembly/SRRxxxxxx/contigs.fasta
#
# Output:
#   work/quast/SRRxxxxxx/
# ============================================================

set -euo pipefail

module purge
module load StdEnv/2020
module load gcc/9.3.0
module load quast/5.0.2

BASE="${KP_WORKDIR:-$PWD/work}"
LIST="${BASE}/srr_list.txt"
ASSEMBLY="${BASE}/assembly"
QUAST_OUT="${BASE}/quast"
STATUS="${BASE}/logs/quast_status"

mkdir -p "$QUAST_OUT" "$STATUS"

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [[ -z "${SRR}" ]]; then
    echo "ERROR: No accession found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

CONTIG="${ASSEMBLY}/${SRR}/contigs.fasta"
OUT="${QUAST_OUT}/${SRR}"
REPORT="${OUT}/report.tsv"

echo "============================================================"
echo "Sample: ${SRR}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "============================================================"

# Skip completed QUAST analyses.
if [[ -s "$REPORT" ]]; then
    echo "${SRR}: QUAST report already exists; skipping."
    rm -f "${STATUS}/${SRR}.failed"
    touch "${STATUS}/${SRR}.complete"
    exit 0
fi

if [[ ! -s "$CONTIG" ]]; then
    echo "ERROR: Assembly not found for ${SRR}"
    echo "Expected: $CONTIG"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

rm -rf "$OUT"

quast.py \
    "$CONTIG" \
    --output-dir "$OUT" \
    --threads "${SLURM_CPUS_PER_TASK}"

if [[ ! -s "$REPORT" ]]; then
    echo "ERROR: QUAST report was not created for ${SRR}"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

rm -f "${STATUS}/${SRR}.failed"
touch "${STATUS}/${SRR}.complete"

echo "Completed QUAST successfully for ${SRR}"
echo "Report: $REPORT"
echo "Finished: $(date)"

