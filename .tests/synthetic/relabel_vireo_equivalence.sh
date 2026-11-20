#!/usr/bin/env bash
# Old-vs-new equivalence harness for workflow/scripts/relabel_vireo.R.
#
# Runs relabel_vireo.R from two git refs on synthetic vireo outputs (from
# make_vireo_outputs.py) inside the pipeline's R container, with a mocked
# Snakemake S4 object, then compares every output: plain files byte-for-byte,
# gzipped ones after decompression. Also records each run's exit status, so a
# scenario that must fail (missing2) is checked to fail under both refs.
#
# Usage: relabel_vireo_equivalence.sh OLD_REF NEW_REF [WORKDIR]
#   (run from inside the snp_demux repo; needs singularity + the cached r image)
set -euo pipefail
old=${1:?old ref}; new=${2:?new ref}
repo=$(git rev-parse --show-toplevel)
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
W=${3:-$(mktemp -d)}; rm -rf "$W"; mkdir -p "$W"
IMG=${R_IMG:-/g/data/a56/software/singularity/images/$(printf %s docker://ghcr.io/swarbricklab/snp_demux-r:r4.2 | md5sum | cut -d' ' -f1).simg}

python3 "$here/make_vireo_outputs.py" "$W/data" --cells "${CELLS:-25000}"
git -C "$repo" show "$old:workflow/scripts/relabel_vireo.R" > "$W/old.R"
git -C "$repo" show "$new:workflow/scripts/relabel_vireo.R" > "$W/new.R"

# Mocked `snakemake` object: just the slots relabel_vireo.R reads.
cat > "$W/run.R" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)  # script, scenario_dir, out_dir
d <- args[2]; o <- args[3]; dir.create(o, showWarnings = FALSE, recursive = TRUE)
# A real one-page PDF and an empty placeholder, as the vireo rule leaves them.
pdf(file.path(d, "est.pdf")); plot(1:10); invisible(dev.off())
file.create(file.path(d, "input.pdf"))
setClass("Snakemake", representation(input = "list", output = "list", wildcards = "list", log = "list"))
snakemake <- new("Snakemake",
  input = list(samples_csv = file.path(d, "samples.csv"), donor_ids = file.path(d, "donor_ids.tsv"),
               prob_doublet = file.path(d, "prob_doublet.tsv.gz"), prob_singlet = file.path(d, "prob_singlet.tsv.gz"),
               prop_ambient = file.path(d, "prop_ambient.tsv"), fig_est = file.path(d, "est.pdf"),
               fig_input = file.path(d, "input.pdf"), summary = file.path(d, "summary.tsv")),
  output = list(donor_ids = file.path(o, "donor_ids.tsv"), sample_map = file.path(o, "sample_map.tsv"),
                prob_doublet = file.path(o, "prob_doublet.tsv.gz"), prob_singlet = file.path(o, "prob_singlet.tsv.gz"),
                prop_ambient = file.path(o, "prop_ambient.tsv"), fig_est_png = file.path(o, "fig_est.png"),
                input_dist_png = file.path(o, "fig_input.png"), summary = file.path(o, "summary.tsv")),
  wildcards = list(capture_id = "synth_cap"),
  log = list(file.path(o, "relabel.log")))
t0 <- Sys.time()
source(args[1])
cat(sprintf("elapsed %.1f s\n", as.numeric(Sys.time() - t0, units = "secs")), file = file.path(o, "elapsed.txt"))
EOF

status=0
for sc in missing1 all_known zero_cells missing2; do
    for v in old new; do
        rc=0
        singularity exec --bind "$W" "$IMG" Rscript "$W/run.R" "$W/$v.R" "$W/data/$sc" "$W/out/$v/$sc" \
            > "$W/out.$v.$sc.txt" 2>&1 || rc=$?
        echo "$rc" > "$W/rc.$v.$sc"
    done
    rc_old=$(cat "$W/rc.old.$sc"); rc_new=$(cat "$W/rc.new.$sc")
    echo "== $sc: exit old=$rc_old new=$rc_new  ($(cat "$W/out/old/$sc/elapsed.txt" 2>/dev/null || echo '-') -> $(cat "$W/out/new/$sc/elapsed.txt" 2>/dev/null || echo '-'))"
    [[ $rc_old == "$rc_new" ]] || { echo "   FAIL exit status differs"; status=1; }
    for f in donor_ids.tsv sample_map.tsv prob_doublet.tsv.gz prob_singlet.tsv.gz prop_ambient.tsv summary.tsv fig_est.png fig_input.png; do
        a=$W/out/old/$sc/$f; b=$W/out/new/$sc/$f
        [[ -e $a || -e $b ]] || continue
        if [[ $f == *.gz ]]; then cmp -s <(zcat "$a") <(zcat "$b"); else cmp -s "$a" "$b"; fi \
            && echo "   identical  $f" || { echo "   DIFFERENT  $f"; status=1; }
    done
    # Flag NA donor labels in any relabelled header (pre-existing edge case check)
    for f in prob_singlet.tsv.gz prob_doublet.tsv.gz prop_ambient.tsv; do
        b=$W/out/new/$sc/$f; [[ -e $b ]] || continue
        zcat -f "$b" | head -1 | tr '\t,' '\n\n' | grep -qx NA && echo "   NOTE  'NA' donor label in $f header (new; old: $(zcat -f "$W/out/old/$sc/$f" | head -1 | tr '\t,' '\n\n' | grep -cx NA) NA)"
    done
    [[ $rc_new != 0 ]] && echo "   error: $(grep -h -m1 -o 'More than one missing donor[^"]*\|Error[^"]*' "$W/out/new/$sc/relabel.log" "$W/out.new.$sc.txt" 2>/dev/null | head -1)"
done
echo "workdir: $W"
exit $status
