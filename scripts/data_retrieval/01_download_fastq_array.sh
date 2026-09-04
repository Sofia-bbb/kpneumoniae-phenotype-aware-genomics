#!/bin/bash
#SBATCH --job-name=KP_download
#SBATCH --account=def-walshj
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/download_%A_%a.out
#SBATCH --error=logs/download_%A_%a.err

set -euo pipefail

module load sra-toolkit/3.0.9

ACCESSION_FILE="${KP_ACCESSION_FILE:-$PWD/config/sra_accessions_1181.txt}"
WORKDIR="${KP_WORKDIR:-$PWD/work}"\nSRA_DIR="${WORKDIR}/sra_cache"
FASTQ_DIR="${WORKDIR}/raw_fastq"
TMP_BASE="${WORKDIR}/tmp"

ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCESSION_FILE")
TASK_TMP="${TMP_BASE}/${ACC}_${SLURM_JOB_ID}"

if [[ -z "$ACC" ]]; then
    echo "ERROR: No accession found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$TASK_TMP"

echo "Starting accession: $ACC"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"

# Skip accessions already downloaded and compressed successfully
if [[ -s "${FASTQ_DIR}/${ACC}_1.fastq.gz" && \
      -s "${FASTQ_DIR}/${ACC}_2.fastq.gz" ]]; then
    echo "$ACC is already complete. Skipping."
    rm -rf "$TASK_TMP"
    exit 0
fi

# Download only when a valid cached SRA file is not already available
SRA_FILE="${SRA_DIR}/${ACC}/${ACC}.sra"

if [[ -s "$SRA_FILE" ]]; then
    echo "Cached SRA file found for $ACC. Validating it."
    vdb-validate "$SRA_FILE"
else
    echo "No cached SRA file found. Downloading $ACC."
    prefetch "$ACC" \
        --output-directory "$SRA_DIR" \
        --max-size 20G

    SRA_FILE="${SRA_DIR}/${ACC}/${ACC}.sra"

    if [[ ! -s "$SRA_FILE" ]]; then
        echo "ERROR: SRA download was not created for $ACC"
        exit 1
    fi

    vdb-validate "$SRA_FILE"
fi

# Convert to paired-end FASTQ
fasterq-dump "$SRA_FILE" \
    --split-files \
    --outdir "$FASTQ_DIR" \
    --temp "$TASK_TMP" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --progress

# Confirm that both paired FASTQ files were created
if [[ ! -s "${FASTQ_DIR}/${ACC}_1.fastq" || \
      ! -s "${FASTQ_DIR}/${ACC}_2.fastq" ]]; then
    echo "ERROR: Paired FASTQ files were not created for $ACC"
    exit 1
fi

# Compress the FASTQ files
gzip -f "${FASTQ_DIR}/${ACC}_1.fastq"
gzip -f "${FASTQ_DIR}/${ACC}_2.fastq"

# Validate the compressed files
gzip -t "${FASTQ_DIR}/${ACC}_1.fastq.gz"
gzip -t "${FASTQ_DIR}/${ACC}_2.fastq.gz"

# Remove temporary and cached SRA data only after successful validation
rm -rf "$TASK_TMP"
rm -rf "${SRA_DIR:?}/${ACC}"

echo "Successfully completed: $ACC"
echo "Date: $(date)"

