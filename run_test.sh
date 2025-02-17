#! /bin/bash

set -e
eval  "$(conda shell.bash hook)"
conda activate snakemake_7.32.4

if [[ "$(hostname)" == *"nci"* ]]; then
    global_profile="--profile profiles/global/nci_a56"
    workflow_profile="--workflow-profile profiles/workflow "
    module load singularity
    mkdir -p logs/joblogs
else
    global_profile=""
fi

snakemake $global_profile $workflow_profile \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml \
    $@
