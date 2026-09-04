#!/bin/bash
#SBATCH --job-name=KP_Trimmomatic
#SBATCH --account=def-walshj
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-300%20
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# ============================================================
# Script: 03_Trimmomatic.sh
# Project: Klebsiella pneumoniae genotype-phenotype study
# Authors: Sofia B. Mohamed and Jacob P. Walsh
#
# Purpose:
#   Trim adapters and low-quality bases from paired-end
#   Illumina FASTQ files using Trimmomatic.
#
# Input:
#   work/raw_fastq/
#       SRRxxxxxx_1.fastq.gz
#       SRRxxxxxx_2.fastq.gz
#
# Output:
#   work/trimmed_fastq/
#       SRRxxxxxx_1_paired.fastq.gz
#       SRRxxxxxx_2_paired.fastq.gz
#       SRRxxxxxx_1_unpaired.fastq.gz
#       SRRxxxxxx_2_unpaired.fastq.gz
#
# Trimming parameters:
#   ILLUMINACLIP: TruSeq3 paired-end adapter removal
#   LEADING:3
#   TRAILING:3
#   SLIDINGWINDOW:4:15
#   MINLEN:36
# ============================================================

set -euo pipefail

module load StdEnv/2023
module load trimmomatic

BASE="${KP_WORKDIR:-$PWD/work}"
LIST="${BASE}/srr_list.txt"

RAW="${BASE}/raw_fastq"
TRIM="${BASE}/trimmed_fastq"
STATUS="${BASE}/logs/trimmomatic_status"
TRIMLOG="${BASE}/logs/trimmomatic_reports"

mkdir -p "$RAW" "$TRIM" "$STATUS" "$TRIMLOG"

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [[ -z "${SRR}" ]]; then
    echo "ERROR: No accession found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

R1="${RAW}/${SRR}_1.fastq.gz"
R2="${RAW}/${SRR}_2.fastq.gz"

P1="${TRIM}/${SRR}_1_paired.fastq.gz"
U1="${TRIM}/${SRR}_1_unpaired.fastq.gz"
P2="${TRIM}/${SRR}_2_paired.fastq.gz"
U2="${TRIM}/${SRR}_2_unpaired.fastq.gz"

LOG="${TRIMLOG}/${SRR}_trimmomatic.log"
SUMMARY="${TRIMLOG}/${SRR}_summary.txt"

echo "============================================================"
echo "Sample: ${SRR}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "============================================================"

# Skip samples already completed successfully.
if [[ -s "$P1" && -s "$P2" ]] &&
   gzip -t "$P1" "$P2" 2>/dev/null; then
    echo "${SRR}: Valid paired trimmed FASTQ files already exist; skipping."
    rm -f "${STATUS}/${SRR}.failed"
    touch "${STATUS}/${SRR}.complete"
    exit 0
fi

# Confirm both raw paired FASTQ files exist.
if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    echo "ERROR: Missing raw paired FASTQ file(s) for ${SRR}"
    echo "Expected:"
    echo "  $R1"
    echo "  $R2"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

# Validate compressed input files.
gzip -t "$R1"
gzip -t "$R2"

# Locate the TruSeq3 paired-end adapter file supplied with Trimmomatic.
ADAPTERS=$(find "${EBROOTTRIMMOMATIC}" \
    -type f \
    -name "TruSeq3-PE*.fa" \
    2>/dev/null | head -n 1)

if [[ -z "$ADAPTERS" || ! -s "$ADAPTERS" ]]; then
    echo "ERROR: TruSeq3 paired-end adapter file was not found."
    echo "EBROOTTRIMMOMATIC=${EBROOTTRIMMOMATIC:-not_set}"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

echo "Adapter file: $ADAPTERS"

# Remove incomplete outputs from an earlier failed attempt.
rm -f "$P1" "$U1" "$P2" "$U2" "$LOG" "$SUMMARY"

# Run paired-end trimming.
java -jar "$EBROOTTRIMMOMATIC/trimmomatic-0.41.jar" PE \
    -threads "${SLURM_CPUS_PER_TASK}" \
    -phred33 \
    -trimlog "$LOG" \
    -summary "$SUMMARY" \
    "$R1" "$R2" \
    "$P1" "$U1" \
    "$P2" "$U2" \
    "ILLUMINACLIP:${ADAPTERS}:2:30:10" \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36

# Confirm paired outputs were created.
if [[ ! -s "$P1" || ! -s "$P2" ]]; then
    echo "ERROR: Paired trimmed FASTQ files were not created for ${SRR}"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

# Validate all non-empty compressed outputs.
for FILE in "$P1" "$P2" "$U1" "$U2"; do
    if [[ -s "$FILE" ]]; then
        gzip -t "$FILE"
    fi
done

rm -f "${STATUS}/${SRR}.failed"
touch "${STATUS}/${SRR}.complete"

echo "Completed Trimmomatic successfully for ${SRR}"
echo "Paired R1: $P1"
echo "Paired R2: $P2"
echo "Summary:   $SUMMARY"
echo "Finished: $(date)"




