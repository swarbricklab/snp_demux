# ==================================
# Doublet-rate inputs for downstream doublet detection
# ----------------------------------
# Per-barcode doublet posteriors (consensus call + vireo posteriors), consumed by
# expression-based doublet detection downstream (e.g. scDblFinder in mkobj).
# ==================================


rule doublet_posteriors:
    input:
        samples=config["deps"]["samples"],
        cell_assignment=rules.wrangle_demuxafy_assignments.output.cell_assignment,
        vireo=rules.relabel_vireo.output.donor_ids,
    output:
        csv=out_dir / "assignment/{capture_id}/doublet_posteriors.csv",
    log:
        logs / "assignment/{capture_id}/doublet_posteriors.log",
    container:
        config["containers"]["vireosnp"]
    script:
        "../scripts/doublet_posteriors.py"


# One-row per-capture summary: cross-donor doublet fraction (hard and posterior),
# every declared sample's share of singlets (sum_p_sq), and the expected total /
# within-donor doublet rates derived from them.
rule doublet_summary:
    input:
        samples=config["deps"]["samples"],
        posteriors=rules.doublet_posteriors.output.csv,
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
