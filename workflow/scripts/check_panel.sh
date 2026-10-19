#!/bin/bash
# Preflight check for a per-pool informative donor-genotype panel, enforcing the
# three input requirements in docs/data_prep_requirements.md so a misprepared
# panel is caught up front with a clear message instead of a cryptic popscle /
# demuxlet FATAL mid-run. Read-only; writes a short report and exits non-zero on
# a hard failure.
#
# Invoked via Snakemake's `script:` directive, which injects the rule I/O as bash
# associative arrays (${snakemake_input[...]}, ${snakemake_output[...]}). The `:-`
# fallbacks to positional args also let it be run standalone:
#   check_panel.sh <panel.vcf.gz> <report.txt>
set -euo pipefail
VCF="${snakemake_input[vcf]:-${1:-}}"
REPORT="${snakemake_output[report]:-${2:-}}"
: > "$REPORT"
log() { echo "$@" | tee -a "$REPORT"; }
fail=0

log "# Panel preflight: $VCF"

# --- 1. every autosome chr1-22 carries at least one site --------------------
# popscle dsc-pileup fatals on an autosome whose (split) VCF has zero variants.
present=$(bcftools index -s "$VCF" 2>/dev/null | awk '$3>0{c=$1; sub(/^chr/,"",c); print c}' | sort -u)
missing=""
for n in $(seq 1 22); do
  echo "$present" | grep -qx "$n" || missing="$missing chr$n"
done
if [ -n "$missing" ]; then
  log "FAIL  [autosomes] no sites on:$missing"
  log "      popscle dsc-pileup fatals on a zero-variant autosome; keep >=1 site per chr1-22."
  fail=1
else
  log "PASS  [autosomes] all of chr1-22 carry >=1 site"
fi

# --- 3. every site has an R2 INFO value -------------------------------------
# demuxlet --geno-error-coeff fatals on the first site lacking R2.
total=$(bcftools view -H "$VCF" | wc -l)
no_r2=$(bcftools view -H "$VCF" | awk -F'\t' '$8 !~ /(^|;)R2=/{c++} END{print c+0}')
if [ "$no_r2" -gt 0 ]; then
  log "FAIL  [R2] $no_r2 of $total sites lack an R2 INFO value"
  log "      demuxlet --geno-error-coeff fatals on the first site without R2; fill R2 (e.g. =1 for TYPED sites)."
  fail=1
else
  log "PASS  [R2] all $total sites carry an R2 INFO value"
fi

# --- 2. contig ordering (heuristic; authoritative check is vs the BAM @SQ) ---
order=$(bcftools view -h "$VCF" | awk -F'[=,]' '/^##contig/{for(i=1;i<=NF;i++)if($i=="ID")print $(i+1)}')
auto=$(echo "$order" | grep -E '^chr?([1-9]|1[0-9]|2[0-2])$' || true)
lex=$(echo "$auto" | LC_ALL=C sort)
if [ "$auto" = "$lex" ]; then
  log "PASS  [contig-order] panel ##contig autosomes are in lexicographic order (matches Cell Ranger GRCh38 @SQ)"
else
  log "WARN  [contig-order] panel ##contig order is NOT lexicographic (looks karyotypic)."
  log "      popscle requires VCF ##contig order == BAM @SQ order (Cell Ranger is lexicographic:"
  log "      chr1,chr10,...,chr19,chr2,...). Reheader+sort the panel if the pileup BAM is Cell Ranger."
fi

log ""
if [ "$fail" -ne 0 ]; then
  log "RESULT: FAIL — fix the panel before running (see docs/data_prep_requirements.md)."
  exit 1
fi
log "RESULT: OK"
