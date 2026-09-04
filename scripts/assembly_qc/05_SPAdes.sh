#!/bin/bash
#SBATCH --job-name=KP_SPAdes
#SBATCH --account=def-walshj
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-5%2
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# ============================================================
# Script: 05_SPAdes.sh
# Project: Klebsiella pneumoniae genotype-phenotype study
# Authors: Sofia B. Mohamed and Jacob P. Walsh
#
# Purpose:
#   Assemble paired-end, quality-trimmed Illumina reads using
#   SPAdes in isolate mode.
#
# Input:
#   work/trimmed_fastq/
#       SRRxxxxxx_1_paired.fastq.gz
#       SRRxxxxxx_2_paired.fastq.gz
#
# Output:
#   work/assembly/SRRxxxxxx/
#       contigs.fasta
#       scaffolds.fasta
#       spades.log
#
# Notes:
#   One SLURM task processes one BioSample representative run.
#   Only paired reads are used for consistent assembly.
# ============================================================

set -euo pipefail

module load StdEnv/2023
module load spades

BASE="${KP_WORKDIR:-$PWD/work}"
LIST="${BASE}/srr_list.txt"

TRIM="${BASE}/trimmed_fastq"
ASSEMBLY="${BASE}/assembly"
STATUS="${BASE}/logs/spades_status"

mkdir -p "$ASSEMBLY" "$STATUS"

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [[ -z "${SRR}" ]]; then
    echo "ERROR: No accession found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

R1="${TRIM}/${SRR}_1_paired.fastq.gz"
R2="${TRIM}/${SRR}_2_paired.fastq.gz"

OUT="${ASSEMBLY}/${SRR}"
CONTIGS="${OUT}/contigs.fasta"
SCAFFOLDS="${OUT}/scaffolds.fasta"

echo "============================================================"
echo "Sample: ${SRR}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Threads: ${SLURM_CPUS_PER_TASK}"
echo "Memory limit: 32 GB"
echo "============================================================"

# Skip a sample that already has a completed, non-empty assembly.
if [[ -s "$CONTIGS" ]]; then
    echo "${SRR}: Existing contigs.fasta found; skipping."
    rm -f "${STATUS}/${SRR}.failed"
    touch "${STATUS}/${SRR}.complete"
    exit 0
fi

# Confirm that both trimmed paired-read files exist.
if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    echo "ERROR: Missing trimmed paired FASTQ file(s) for ${SRR}"
    echo "Expected:"
    echo "  $R1"
    echo "  $R2"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

# Validate the compressed input files.
gzip -t "$R1"
gzip -t "$R2"

# If a valid interrupted SPAdes run exists, attempt to continue it.
if [[ -d "$OUT" && -s "${OUT}/params.txt" ]]; then

    echo "Interrupted SPAdes directory detected; attempting --continue."

    spades.py \
        --continue \
        -o "$OUT"

else

    # Remove an incomplete output directory that cannot be resumed.
    rm -rf "$OUT"

    spades.py \
        --isolate \
        -1 "$R1" \
        -2 "$R2" \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --memory 32 \
        -o "$OUT"
fi

# Confirm that the principal assembly output exists.
if [[ ! -s "$CONTIGS" ]]; then
    echo "ERROR: SPAdes did not produce contigs.fasta for ${SRR}"
    touch "${STATUS}/${SRR}.failed"
    exit 1
fi

# Record basic assembly information.
CONTIG_COUNT=$(grep -c '^>' "$CONTIGS")
ASSEMBLY_SIZE=$(awk '
    /^>/ {next}
    {total += length($0)}
    END {print total}
' "$CONTIGS")

echo "Contigs produced: ${CONTIG_COUNT}"
echo "Assembly size: ${ASSEMBLY_SIZE} bp"

rm -f "${STATUS}/${SRR}.failed"
touch "${STATUS}/${SRR}.complete"

echo "Completed SPAdes successfully for ${SRR}"
echo "Contigs: ${CONTIGS}"

if [[ -s "$SCAFFOLDS" ]]; then
    echo "Scaffolds: ${SCAFFOLDS}"
fi

echo "Finished: $(date)"

