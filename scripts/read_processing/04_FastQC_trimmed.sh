#!/bin/bash
#SBATCH --job-name=KP_FastQC_trimmed
#SBATCH --account=def-walshj
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --array=1-5%2
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# ============================================================
# Script: 04_FastQC_trimmed.sh
# Project: Klebsiella pneumoniae genotype-phenotype study
# Authors: Sofia B. Mohamed and Jacob P. Walsh
#
# Purpose:
#   Perform FastQC quality assessment on paired-end
#   trimmed FASTQ files.
#
# Input:
#   work/trimmed_fastq/
#
# Output:
#   work/trimmed_fastqc/
# ============================================================

set -euo pipefail

module load StdEnv/2023
module load fastqc

BASE="${KP_WORKDIR:-$PWD/work}"
LIST="${BASE}/srr_list.txt"

TRIM="${BASE}/trimmed_fastq"
OUT="${BASE}/trimmed_fastqc"
STATUS="${BASE}/logs/fastqc_trimmed_status"

mkdir -p "$TRIM" "$OUT" "$STATUS"

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [[ -z "$SRR" ]]; then
    echo "No accession found."
    exit 1
fi

R1="${TRIM}/${SRR}_1_paired.fastq.gz"
R2="${TRIM}/${SRR}_2_paired.fastq.gz"

ZIP1="${OUT}/${SRR}_1_paired_fastqc.zip"
ZIP2="${OUT}/${SRR}_2_paired_fastqc.zip"

HTML1="${OUT}/${SRR}_1_paired_fastqc.html"
HTML2="${OUT}/${SRR}_2_paired_fastqc.html"

echo "Processing ${SRR}"

if [[ -s "$ZIP1" && -s "$ZIP2" ]]; then
    echo "Already completed."
    touch "${STATUS}/${SRR}.complete"
    exit 0
fi

gzip -t "$R1"
gzip -t "$R2"

fastqc \
    --threads ${SLURM_CPUS_PER_TASK} \
    --outdir "$OUT" \
    "$R1" "$R2"

if [[ ! -s "$ZIP1" || ! -s "$ZIP2" ]]; then
    echo "FastQC failed."
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

rm -f "${STATUS}/${SRR}.failed"
touch "${STATUS}/${SRR}.complete"

echo "Completed ${SRR}"

