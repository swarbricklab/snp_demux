# ==================================
# Doublet-rate inputs for downstream doublet detection
# ----------------------------------
# Per-capture doublet summary from the consensus calls and the vireo posteriors
# in sample_barcode_assignments.tsv (see add_posteriors.py), consumed by
# expression-based doublet detection downstream (e.g. scDblFinder in mkobj).
# ==================================


# One-row per-capture summary: cross-donor doublet rate (hard and posterior),
# every declared sample's share of singlets (same_donor_pair_prob), and the
# expected total / within-donor doublet rates derived from them.
rule doublet_summary:
    input:
        samples=config["deps"]["samples"],
        cell_assignment=rules.wrangle_demuxafy_assignments.output.cell_assignment,
        assignments=rules.wrangle_demuxafy_assignments.output.sample_barcode_assignments,
    output:
        csv=out_dir / "assignment/{capture_id}/doublet_summary.csv",
    log:
        logs / "assignment/{capture_id}/doublet_summary.log",
    container:
        config["containers"]["vireosnp"]
    params:
        rate_per_1k=get_doublet_rate_per_1k,
    script:
        "../scripts/doublet_summary.py"


rule combine_doublet_summary:
    input:
        summaries=expand(
            out_dir / "assignment/{capture_id}/doublet_summary.csv",
            capture_id=captures.index.tolist(),
        ),
    output:
        csv=out_dir / "assignment/report/doublet_summary.csv",
    log:
        logs / "assignment/report/doublet_summary.log",
    container:
        config["containers"]["vireosnp"]
    script:
        "../scripts/doublet_summary.py"
