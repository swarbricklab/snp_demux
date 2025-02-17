# ==================================
# Dan Roden
# ----------------------------------
# CREATED: 22/09/2020
# MODIFIED: 22/09/2020
# ----------------------------------
# Rules for running Souporcell sample demux
# ------------------------------------------------------------------------------
# ==================================

# ==================================
# Generate cell barcode lists (from )
# ----------------------------------
# TODO: Description - Extract cell subsets based on a meta-data column from infercnv
# NOTE" currently runs one job per SAMPLE not capture, as this streamlines snakemake rules and removes need for checkpoints at this stage
# ==================================
rule pdf_to_png:
    input:
        Singlets_upset_pdf=out_dir/"assignment/{capture_id}/combined_resultsSinglets_upset.pdf"
    output:
        Singlets_upset_png=out_dir/"assignment/{capture_id}/combined_resultsSinglets_upset.png"
    conda:
        "../envs/pdftools.yaml"
    log:
        logs/"pdf_to_png/{capture_id}.log"
    script:
        "../scripts/pdf_to_png.R"
