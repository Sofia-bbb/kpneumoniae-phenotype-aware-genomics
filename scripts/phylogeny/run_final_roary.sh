#!/bin/bash
#SBATCH --job-name=KP_Final_Roary
#SBATCH --account=def-walshj
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/Roary_%j.out
#SBATCH --error=logs/Roary_%j.err

module purge
module load StdEnv/2020
module load gcc/9.3.0
module load roary/3.13.0

cd "${KP_WORKDIR:-$PWD/work}"

mkdir -p "${KP_ROARY_OUTDIR:-phylogeny/roary_output}"

roary \
  -e \
  -n \
  -v \
  -p "${SLURM_CPUS_PER_TASK:-32}" \
  --group_limit 100000 \
  -f "${KP_ROARY_OUTDIR:-phylogeny/roary_output}" \
  15_Final_Integrated_Phylogeny/GFF/*.gff



