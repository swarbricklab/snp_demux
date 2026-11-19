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
