#! /bin/bash

set -e
eval  "$(conda shell.bash hook)"
conda activate snakemake_7.32.4

if [[ "$(hostname)" == *"nci"* || "$(hostname)" == *"gadi"* ]]; then
    # The global profile reads $PROJECT (default to a56 for this lab's storage).
    export PROJECT="${PROJECT:-a56}"
    global_profile="--profile profiles/global/nci"
    workflow_profile="--workflow-profile profiles/workflow "
    module load singularity
    mkdir -p logs/joblogs
else
    global_profile=""
fi

snakemake $global_profile $workflow_profile \
    --snakefile workflow/Snakefile \
    --configfile config/test_config.yaml \
    $@
