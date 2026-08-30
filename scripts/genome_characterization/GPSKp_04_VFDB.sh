#!/bin/bash
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --job-name=GPSKp_VFDB
#SBATCH --output=logs/GPSKp_VFDB_%A_%a.out
#SBATCH --error=logs/GPSKp_VFDB_%A_%a.err

set -euo pipefail

module purge
module load apptainer/1.4.5

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

MANIFEST=$BASE/data/assemblies.tsv
DB=$BASE/data/vfdb
CONTAINER=$BASE/containers/kleborate.sif

LINE=$((SLURM_ARRAY_TASK_ID+1))

SRR=$(sed -n "${LINE}p" "$MANIFEST" | cut -f1)
ASSEMBLY=$(sed -n "${LINE}p" "$MANIFEST" | cut -f2)

OUTDIR=$BASE/results/vfdb/$SRR
STATUS=$BASE/results/vfdb/status

mkdir -p "$OUTDIR" "$STATUS" "$BASE/logs"

echo "================================="
echo "GPS-Kp VFDB"
echo "Sample: $SRR"
echo "Assembly: $ASSEMBLY"
echo "Started: $(date)"
echo "================================="

apptainer exec "$CONTAINER" \
blastn \
-query "$ASSEMBLY" \
-db "$DB" \
-evalue 1e-10 \
-max_target_seqs 5 \
-num_threads 4 \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
> "$OUTDIR/${SRR}_vfdb_raw.tsv"

touch "$STATUS/$SRR.complete"

echo "Finished: $(date)"


