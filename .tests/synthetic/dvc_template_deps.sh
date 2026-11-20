#!/usr/bin/env bash
# Synthetic consumer repo for the DVC stage template (.dvc/template.yaml).
# Renders the template into a throwaway DVC repo with placeholder dep dirs, commits
# the stage, then perturbs each candidate input and reports whether `dvc status`
# flags the stage. Used to verify #56 (depend on pileups + panels, not BAMs).
#
# Usage: dvc_template_deps.sh <git-ref> [workdir]   (needs dvc on PATH, e.g. env dt3)
set -euo pipefail
ref=${1:?git ref}; repo=$(git rev-parse --show-toplevel)
T=${2:-$(mktemp -d)}; rm -rf "$T"; mkdir -p "$T"; cd "$T"
git init -q && dvc init -q
mkdir -p modules/snp_demux/workflow config/snp_demux data/snp/pileup/cellsnp \
         data/snp/vcf/panel/cap1 data/snp/imputation data/sc/cellranger data/snp/demux logs/snp/demux
echo x > modules/snp_demux/run_mod.sh; echo x > modules/snp_demux/workflow/Snakefile
echo x > config/snp_demux/config.yaml; echo p > data/snp/pileup/cellsnp/a
echo v > data/snp/vcf/panel/cap1/informative.vcf.gz; echo i > data/snp/imputation/i.vcf.gz
echo b > data/sc/cellranger/possorted_genome_bam.bam
(echo "stages:"; git -C "$repo" show "$ref:.dvc/template.yaml") > dvc.yaml
dvc commit -f -q
check() { # $1 = file to perturb, $2 = expectation (stale|fresh)
    echo "$RANDOM" >> "$1"
    if dvc status -q snp_demux; then got=fresh; else got=stale; fi
    printf '%-50s expect=%-5s got=%-5s %s\n' "$1" "$2" "$got" "$([[ $got == "$2" ]] && echo PASS || echo FAIL)"
    dvc commit -f -q
}
check data/sc/cellranger/possorted_genome_bam.bam fresh
check data/snp/pileup/cellsnp/a stale
check data/snp/vcf/panel/cap1/informative.vcf.gz stale
