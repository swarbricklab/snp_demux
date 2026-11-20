#!/usr/bin/env bash
# Run workflow/scripts/wrangle_assignments.R on a consensus calls table, outside
# Snakemake (mocked S4 snakemake object), in the pipeline's R container.
#   wrangle_calls.sh CALLS.tsv SAMPLES.csv CAPTURE_ID OUTDIR
# Writes OUTDIR/sample_barcode_assignments.tsv and OUTDIR/cell_assignment.tsv.
set -euo pipefail
calls=$(realpath "$1"); samples=$(realpath "$2"); cap=$3; out=$(realpath -m "$4")
repo=$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)
IMG=${R_IMG:-/g/data/a56/software/singularity/images/$(printf %s docker://ghcr.io/swarbricklab/snp_demux-r:r4.2 | md5sum | cut -d' ' -f1).simg}
mkdir -p "$out"
cat > "$out/run.R" <<RS
setClass("Snakemake", representation(input = "list", output = "list", params = "list", wildcards = "list", log = "list"))
snakemake <- new("Snakemake",
  input = list(sample_map = "$samples", demuxafy_assignments = "$calls"),
  output = list(sample_barcode_assignments = "$out/sample_barcode_assignments.tsv",
                cell_assignment = "$out/cell_assignment.tsv"),
  params = list(assignment_column = "AtLeastHalfSinglet_Individual_Assignment",
                droplettype_column = "AtLeastHalfSinglet_DropletType"),
  wildcards = list(capture_id = "$cap"), log = list("$out/wrangle.log"))
source("$repo/workflow/scripts/wrangle_assignments.R")
RS
singularity exec --bind "$repo,$out,$(dirname "$calls"),$(dirname "$samples")" "$IMG" Rscript "$out/run.R"
