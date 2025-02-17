#! /bin/bash

set -e
eval  "$(conda shell.bash hook)"
conda activate snakemake_7.32.4

if [[ "$(hostname)" == *"nci"* ]]; then
    global_profile="--profile modules/demuxafy/profiles/global/nci_a56"
    workflow_profile="--workflow-profile modules/demuxafy/profiles/workflow "
    module load singularity
    mkdir -p logs/joblogs
else
    global_profile=""
fi

snakemake $global_profile $workflow_profile \
    --snakefile modules/demuxafy/workflow/Snakefile \
    --configfile config/demuxafy/config.yaml \
    $@
