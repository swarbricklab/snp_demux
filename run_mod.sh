#! /bin/bash

set -e
eval  "$(conda shell.bash hook)"
conda activate snakemake_7.32.4

if [[ "$(hostname)" == *"nci"* ]]; then
    global_profile="--profile modules/snp_demux/profiles/global/nci_a56"
    workflow_profile="--workflow-profile modules/snp_demux/profiles/workflow "
    module load singularity
    mkdir -p logs/joblogs
else
    global_profile=""
fi

snakemake $global_profile $workflow_profile \
    --snakefile modules/snp_demux/workflow/Snakefile \
    --configfile config/snp_demux/config.yaml \
    $@
