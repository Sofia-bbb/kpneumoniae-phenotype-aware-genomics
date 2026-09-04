#!/bin/bash
#SBATCH --job-name=KP_FastQC_raw
#SBATCH --account=def-walshj
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --array=1-1188%20
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# ============================================================
# Script: 02_FastQC_raw.sh
# Project: Klebsiella pneumoniae genotype–phenotype study
# Authors: Sofia B. Mohamed and Jacob P. Walsh
# Purpose:
#   Perform FastQC assessment of paired-end raw FASTQ files.
#
# Input:
#   work/raw_fastq/
#       SRRxxxxxx_1.fastq.gz
#       SRRxxxxxx_2.fastq.gz
#
# Output:
#   work/raw_fastqc/
#
# Notes:
#   One SLURM array task processes one representative SRA run.
#   Completed samples are skipped automatically.
# ============================================================

set -euo pipefail

module load StdEnv/2023
module load fastqc

BASE="${KP_WORKDIR:-$PWD/work}"
LIST="${BASE}/srr_list.txt"
RAW="${BASE}/raw_fastq"
OUT="${BASE}/raw_fastqc"
STATUS="${BASE}/logs/fastqc_raw_status"

mkdir -p "$RAW" "$OUT" "$STATUS"

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [[ -z "${SRR}" ]]; then
    echo "ERROR: No accession found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

R1="${RAW}/${SRR}_1.fastq.gz"
R2="${RAW}/${SRR}_2.fastq.gz"

ZIP1="${OUT}/${SRR}_1_fastqc.zip"
ZIP2="${OUT}/${SRR}_2_fastqc.zip"
HTML1="${OUT}/${SRR}_1_fastqc.html"
HTML2="${OUT}/${SRR}_2_fastqc.html"

echo "============================================================"
echo "Sample: ${SRR}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "============================================================"

# Skip samples already processed successfully.
if [[ -s "$ZIP1" && -s "$ZIP2" && -s "$HTML1" && -s "$HTML2" ]]; then
    echo "${SRR}: FastQC outputs already exist; skipping."
    rm -f "${STATUS}/${SRR}.failed"
    touch "${STATUS}/${SRR}.complete"
    exit 0
fi

# Confirm that both paired FASTQ files exist.
if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    echo "ERROR: Missing paired FASTQ file(s) for ${SRR}"
    echo "Expected:"
    echo "  $R1"
    echo "  $R2"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

# Confirm that both compressed FASTQ files are valid.
gzip -t "$R1"
gzip -t "$R2"

# Run FastQC.
fastqc \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --outdir "$OUT" \
    "$R1" "$R2"

# Confirm that all expected FastQC outputs were created.
if [[ ! -s "$ZIP1" || ! -s "$ZIP2" || \
      ! -s "$HTML1" || ! -s "$HTML2" ]]; then
    echo "ERROR: FastQC outputs are incomplete for ${SRR}"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

rm -f "${STATUS}/${SRR}.failed"
touch "${STATUS}/${SRR}.complete"

echo "Completed FastQC successfully for ${SRR}"
echo "Finished: $(date)"

