#!/bin/bash
# Apply the demux-side imputation-quality cut to a per-pool donor panel.
#
# The precision/recall trade-off for demultiplexing is made here, at read time,
# where its effect on the demux calls can be measured. Filtering is
# non-destructive: the upstream panel is untouched, and a different threshold is
# just a re-filter. This relies on upstream imputation being permissive and
# carrying INFO/R2 through, which is swarbricklab/imputation#24 and not yet
# done: today upstream still filters at its own post_r2 (default 0.8), so a
# min_r2 below that recovers nothing.
#
#   * imputed sites with INFO/R2 < min_r2 are dropped;
#   * directly genotyped sites are kept unconditionally ("trust measured sites,
#     add only high-confidence imputed ones");
#   * if min_maf > 0, sites with INFO/MAF < min_maf are dropped as well.
#
# "Imputed" is read from Minimac's flags: INFO/IMPUTED set, or (if the header
# declares only TYPED) INFO/TYPED not set. With neither flag declared every site
# is treated as imputed. A site with no R2 value is kept. A panel without an
# INFO/R2 header is passed through unchanged, with a warning.
#
# Invoked via Snakemake's `script:` (bash associative arrays), or standalone:
#   filter_panel.sh <in.vcf.gz> <out.vcf.gz> <min_r2> <min_maf>
set -euo pipefail
if declare -p snakemake_input >/dev/null 2>&1; then
    IN="${snakemake_input[vcf]}"
    OUT="${snakemake_output[vcf]}"
    MIN_R2="${snakemake_params[min_r2]}"
    MIN_MAF="${snakemake_params[min_maf]}"
else
    IN="$1"
    OUT="$2"
    MIN_R2="${3:-0.8}"
    MIN_MAF="${4:-0}"
fi

header=$(bcftools view -h "$IN")
declares() { grep -q "^##INFO=<ID=$1," <<<"$header"; }
positive() { awk -v x="$1" 'BEGIN { exit !(x > 0) }'; }

terms=()
if positive "$MIN_R2"; then
    if declares R2; then
        if declares IMPUTED; then
            terms+=("(INFO/IMPUTED=1 && INFO/R2<$MIN_R2)")
        elif declares TYPED; then
            terms+=("(INFO/TYPED=0 && INFO/R2<$MIN_R2)")
        else
            terms+=("INFO/R2<$MIN_R2")
        fi
    else
        echo "WARNING: $IN declares no INFO/R2; imputation R2 filter skipped" >&2
    fi
fi
if positive "$MIN_MAF"; then
    if ! declares MAF; then
        echo "ERROR: min_maf=$MIN_MAF but $IN declares no INFO/MAF" >&2
        exit 1
    fi
    terms+=("INFO/MAF<$MIN_MAF")
fi

if [ "${#terms[@]}" -eq 0 ]; then
    expr=""
    bcftools view -Oz -o "$OUT" "$IN"
else
    expr=$(printf ' || %s' "${terms[@]}")
    expr=${expr:4}
    bcftools view -e "$expr" -Oz -o "$OUT" "$IN"
fi
bcftools index -f "$OUT"

n_in=$(bcftools view -H "$IN" | wc -l | tr -d ' ')
n_out=$(bcftools view -H "$OUT" | wc -l | tr -d ' ')
echo "min_r2=$MIN_R2 min_maf=$MIN_MAF exclude='${expr}' sites: $n_in -> $n_out" >&2
