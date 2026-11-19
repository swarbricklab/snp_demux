#!/bin/bash
# Missing-genotype scenarios for the OneK1K fixture (issue #4).
#
# The golden run has a reference genotype for every donor. These scenarios drop
# donors from the per-pool panel to exercise the missing-genotype logic:
#
#   missing1  one donor has no genotype -> vireo learns it, souporcell gives it
#             the leftover cluster; the run must complete and recover that
#             donor's cells (checked against the golden by compare_missing.py)
#   missing2  two donors have no genotype -> unsupported; the run must fail with
#             the explicit "more than one missing donor" error, not something
#             cryptic or, worse, a silent success
#
# Usage (called by run_test.sh --scenario NAME):
#   scenario.sh prepare NAME   derive the panel + config overlay
#   scenario.sh check NAME RC  verify the outcome given snakemake's exit code
#
# Which donors are dropped: $MISSING_DONORS (comma-separated) if set; otherwise
# the donor(s) with the most golden singlets, so the recovered cluster is large
# enough to check; otherwise the first donor(s) in config/samples.csv.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
CAPTURE=onek1k
BCFTOOLS_IMAGE="docker://quay.io/biocontainers/bcftools:1.13--h3a49de5_0"

action="$1"
name="$2"
case "$name" in
missing1) n_missing=1 ;;
missing2) n_missing=2 ;;
*)
    echo "unknown scenario: $name (missing1|missing2)" >&2
    exit 2
    ;;
esac
dir="$HERE/scenarios/$name"

bcftools_() {
    if command -v bcftools >/dev/null 2>&1; then
        bcftools "$@"
    else
        singularity exec --bind "$ROOT" "$BCFTOOLS_IMAGE" bcftools "$@"
    fi
}

pick_donors() {
    if [[ -n "${MISSING_DONORS:-}" ]]; then
        echo "$MISSING_DONORS"
        return
    fi
    local stats="$HERE/golden/assignment/$CAPTURE/all_singlets_stats.tsv"
    if [[ -s "$stats" ]]; then
        # fixture sample_id == donor_id; sort by number_of_cells, largest first
        awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
            {print $h["number_of_cells"]"\t"$h["sample_id"]}' "$stats" |
            sort -k1,1nr | head -n "$n_missing" | cut -f2 | paste -sd, -
    else
        awk -F, -v c="$CAPTURE" 'NR>1 && $3==c{print $2}' "$ROOT/config/samples.csv" |
            head -n "$n_missing" | paste -sd, -
    fi
}

prepare() {
    local src="$ROOT/.tests/integration/data/pileup/panel/$CAPTURE/informative.vcf.gz"
    local donors
    donors=$(pick_donors)
    echo "scenario $name: dropping donor(s) $donors from the panel" >&2
    mkdir -p "$dir/panel/$CAPTURE"
    bcftools_ view -s "^$donors" -Oz -o "$dir/panel/$CAPTURE/informative.vcf.gz" "$src"
    bcftools_ index -f "$dir/panel/$CAPTURE/informative.vcf.gz"
    echo "$donors" >"$dir/missing_donors.txt"
    cat >"$dir/config.yaml" <<YAML
# Overlay on config/test_config.yaml for scenario $name (generated).
deps:
  panel_root: .tests/integration/scenarios/$name/panel
outs:
  root: .tests/integration/scenarios/$name/results
  logs: .tests/integration/scenarios/$name/logs
YAML
    echo "$dir/config.yaml"
}

check() {
    local rc="$1"
    local donors
    donors=$(cat "$dir/missing_donors.txt")
    if [[ "$name" == missing1 ]]; then
        if [[ "$rc" -ne 0 ]]; then
            echo "FAIL: one donor without a genotype ($donors) must be supported, but snakemake exited $rc" >&2
            exit 1
        fi
        python3 "$HERE/compare_missing.py" \
            "$HERE/golden/assignment/$CAPTURE/cell_assignment.tsv" \
            "$dir/results/assignment/$CAPTURE/cell_assignment.tsv" \
            "$donors"
    else
        if [[ "$rc" -eq 0 ]]; then
            echo "FAIL: two donors without a genotype ($donors) is unsupported, but the run succeeded" >&2
            exit 1
        fi
        if grep -rqE "More than one missing donor|UNSUPPORTED number of samples without a known genotype" "$dir/logs"; then
            echo "PASS: run failed with the explicit unsupported-missing-donors error"
        else
            echo "FAIL: run failed, but not with the expected missing-donor error; see $dir/logs" >&2
            exit 1
        fi
    fi
}

case "$action" in
prepare) prepare ;;
check) check "${3:?exit code}" ;;
*)
    echo "usage: scenario.sh prepare|check NAME [RC]" >&2
    exit 2
    ;;
esac
