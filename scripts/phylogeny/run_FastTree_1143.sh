#!/bin/bash
#SBATCH --job-name=KP_FT_1143
#SBATCH --account=def-walshj
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/FastTree_1143_%j.out
#SBATCH --error=logs/FastTree_1143_%j.err

set -euo pipefail

cd "${KP_WORKDIR:-$PWD/work}"

module purge
module load StdEnv/2020
module load fasttree/2.1.11

ALIGNMENT="${KP_CORE_ALIGNMENT:-phylogeny/roary_output/core_gene_alignment.aln}"

OUTDIR="${KP_FASTTREE_OUTDIR:-phylogeny/fasttree}"
mkdir -p "$OUTDIR"
TREE="${OUTDIR}/KP_Final_1143_core_genome_tree.nwk"

echo "===== FINAL FASTTREE ====="
echo "Start: $(date)"
echo "Alignment: $ALIGNMENT"
echo "Expected isolates: 1143"

echo
echo "===== ALIGNMENT SEQUENCE COUNT ====="
grep -c '^>' "$ALIGNMENT"

FastTree \
-nt \
-gtr \
-gamma \
"$ALIGNMENT" \
> "$TREE"

echo
echo "===== TREE CREATED ====="
ls -lh "$TREE"

echo
echo "===== TREE TIP COUNT ====="

python3 - "$TREE" <<'PY'
import re
import sys

tree = open(sys.argv[1]).read().strip()

tips = re.findall(r'(?<=[(,])([^():,;]+)(?=:)', tree)

print("Tree tips:", len(tips))
print("Unique tips:", len(set(tips)))

if len(tips) != 1143:
    raise SystemExit(
        f"ERROR: expected 1143 tips but detected {len(tips)}"
    )
PY

echo
echo "Completed: $(date)"


