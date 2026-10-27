#! /bin/bash
# Run the snp_demux demultiplexing pipeline on the bundled OneK1K test fixture
# (all four methods + the Demuxafy consensus), locally with singularity, and
# (with --check) verify the fresh output against the DVC-tracked golden.
#
# Assumes:
#   * conda env `snakemake_7.32.4`
#   * the fixture + golden are materialised:  dvc pull   (or dvc checkout)
#   * the pipeline's singularity images + genome FASTA exist at the paths in
#     config/test_config.yaml (NCI-resident; cached in the shared prefix)
#
# Usage:
#   ./run_test.sh            # run all methods into .tests/integration/results
#   ./run_test.sh --check    # run, then content-compare .tests/integration/results to .tests/integration/golden
#   ./run_test.sh -n         # dry run
#   CORES=8 ./run_test.sh    # more local cores
set -euo pipefail

# Make conda usable on a bare (non-login) compute shell, e.g. under `qxub exec`
# in CI: the binary may be on PATH but its shell functions are not set up, so
# `conda activate` errors with "Run 'conda init'". Source the profile if conda is
# absent, then always install the shell hook before activating. Override
# CONDA_PROFILE if your install lives elsewhere.
if ! command -v conda >/dev/null 2>&1; then
    source "${CONDA_PROFILE:-/g/data/a56/conda/miniconda/etc/profile.d/conda.sh}"
fi
eval "$(conda shell.bash hook)"
conda activate snakemake_7.32.4

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# The global profile reads $PROJECT; the resource keys (project/storage) are
# harmless no-ops under local execution but keep the config consistent.
export PROJECT="${PROJECT:-a56}"
module load singularity 2>/dev/null || true

# Pull --check out of the args passed through to snakemake.
check=false
sm_args=()
for a in "$@"; do
    if [[ "$a" == "--check" ]]; then check=true; else sm_args+=("$a"); fi
done

# Local execution with singularity (not the PBS fan-out profile): run_test.sh is
# for the small CI fixture and is invoked inside a single compute job in CI, which
# cannot itself submit cluster jobs. Containers are reused from the shared prefix.
snakemake \
    --snakefile "$HERE/workflow/Snakefile" \
    --configfile "$HERE/config/test_config.yaml" \
    --workflow-profile "$HERE/profiles/workflow" \
    --use-singularity \
    --singularity-args "--bind /g/data/a56,/scratch/a56" \
    --singularity-prefix "${SINGULARITY_PREFIX:-/g/data/a56/software/singularity/images}" \
    --cores "${CORES:-4}" \
    ${sm_args[@]+"${sm_args[@]}"}

if [[ "$check" == true ]]; then
    echo "=== comparing .tests/integration/results to the golden (.tests/integration/golden) ==="
    python3 "$HERE/.tests/integration/compare_golden.py" "$HERE/.tests/integration/golden" "$HERE/.tests/integration/results"
fi
