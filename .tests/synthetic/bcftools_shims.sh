#!/usr/bin/env bash
# Create PATH shims for bcftools/bgzip/tabix that run the pipeline's *cached*
# bcftools image (no registry pull, so they work on offline compute nodes).
# Used to (a) run the bcftools-backed unit tests that otherwise skip, and
# (b) run .tests/integration/scenario.sh, whose docker:// fallback needs internet.
#
# Usage: eval "$(.tests/synthetic/bcftools_shims.sh DIR)"   # prints an export PATH=DIR:$PATH line
set -euo pipefail
dir=${1:?shim dir}; mkdir -p "$dir"
prefix=${SINGULARITY_PREFIX:-/g/data/a56/software/singularity/images}
img="$prefix/$(printf %s docker://quay.io/biocontainers/bcftools:1.13--h3a49de5_0 | md5sum | cut -d' ' -f1).simg"
[[ -s $img ]] || { echo "cached image not found: $img" >&2; exit 1; }
for t in bcftools bgzip tabix; do
    printf '#!/bin/bash\nexec singularity exec --bind /g/data/a56,/scratch/a56 %s %s "$@"\n' "$img" "$t" > "$dir/$t"
    chmod +x "$dir/$t"
done
echo "export PATH=$dir:\$PATH"
