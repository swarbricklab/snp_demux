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
rule generate_upset_plots:
    input:
        result_Demuxlet=demuxlet/"{capture_id}/demuxlet.best",
        result_Freemuxlet=freemuxlet/"{capture_id}/freemuxlet.clust1.samples.gz", 
        # TODO: QUESTION: why is this using a "tmp" input cluster file?
        result_Souporcell=souporcell/"{capture_id}/clusters_tmp.tsv",
        #result_Souporcell=souporcell/"{capture_id}/clusters.tsv",
        result_Vireo=vireo/"{capture_id}/donor_ids.tsv"
    output:
        upset_plot_pool=combined/"{capture_id}/per_pool.pdf",
        upset_plot_per_sample="results/demuxafy/combined/{capture_id}/{sample_id}/per_sample.pdf"
    conda:
        # TODO: add env?
        #"../envs/r_v36.yaml"
    log:
        logs/"combine_results.generate_upset_plots/{capture_id}.log"
    script:
        "../scripts/combine_results.generate_upset_plots.R"
