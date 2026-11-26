#!/bin/bash
# The Scorer against vireo's own output on snp_demux main's OneK1K fixture (14 donors, 1,400 cells): both
# loaders, then validate-vireo against the golden vireo folder. Exits non-zero unless every cell is identical
# in vireo's eight columns and the four vireo files written by --vireo-out match the golden ones as text.
#
#   bash tools/tests/validate_onek1k.sh <main checkout> <output folder> [python command]
#
# <main checkout> is a clone of main with .tests/integration/golden and .tests/integration/data pulled
# (dt pull). The python command should be the vireo image's, for exact reproduction, e.g.
#   "singularity exec -B <paths> <snp_demux-vireo:0.5.7 image> /opt/conda/envs/snp_demux-vireo/bin/python"
set -euo pipefail
MAIN=$(cd "$1" && pwd)
OUT=$2
PY=${3:-python}
TOOLS=$(cd "$(dirname "$0")/.." && pwd)/demux_tools
F=$MAIN/.tests/integration
GOLD=$F/golden/vireo/onek1k
mkdir -p "$OUT"
for loader in native vireo; do
  $PY -I -B "$TOOLS" score --name onek1k --loader $loader --cellsnp "$F/data/pileup/cellsnp/short/onek1k" \
      --vcf "$F/data/pileup/panel/onek1k/informative.vcf.gz" --vireo-out "$OUT/onek1k_$loader.vireo" --force \
      --out "$OUT/onek1k_$loader"
  $PY -I -B "$TOOLS" validate-vireo --scored "$OUT/onek1k_$loader" --vireo-dir "$GOLD" \
      --files "$OUT/onek1k_$loader.vireo" --out "$OUT/onek1k_$loader.validation.json"
done
