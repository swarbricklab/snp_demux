# ==================================
# HTML demux reports
# ----------------------------------
# One page per capture plus a searchable index. Output location comes from
# config['outs']['reports'] so a dataset can publish these alongside its other
# reports; see report_dir in the Snakefile.
# ==================================


def get_samples_in_capture(wildcards):
    """
    Return the sample_ids declared for wildcards.capture_id, comma-separated.
    """
    capture = samples.loc[samples['capture_id'] == wildcards.capture_id]
    return ",".join(capture.index.tolist())


def get_pool_for_capture(wildcards):
    """
    Return the pool_id wildcards.capture_id belongs to, or "" if unrecorded.

    pool_id is optional and used only for display and grouping in the reports;
    no other rule reads it. A blank cell in captures.csv reaches us as a pandas
    NaN, which Snakemake serialises into the generated R script as a bare `nan`
    -- an undefined symbol in R, so the job dies with "object 'nan' not found"
    before the script runs. Coerce to a string so a missing optional field
    cannot fail a rule.
    """
    pool = captures.loc[wildcards.capture_id, "pool_id"]
    if pool is None or (isinstance(pool, float) and pd.isna(pool)):
        return ""
    return str(pool)


rule report_capture_metrics:
    """
    Per-capture roll-up, read by both the capture page and the index.

    Separate from report_capture because Snakemake allows an .Rmd script only a
    single output file. Keeping it separate also means the page and the index
    quote the same numbers instead of each deriving their own.
    """
    input:
        assignments=out_dir/"assignment/{capture_id}/sample_barcode_assignments.tsv",
        cell_assignment=out_dir/"assignment/{capture_id}/cell_assignment.tsv",
        ambient_summary=out_dir/"vireo/{capture_id}/ambient_summary.csv",
        souporcell_ambient=out_dir/"souporcell/{capture_id}/ambient_rna.txt"
    output:
        metrics=report_dir/"metrics/{capture_id}.tsv"
    params:
        declared_samples=get_samples_in_capture,
        pool_id=get_pool_for_capture
    container:
        config["containers"]["r"]
    # conda:  # swap for the container above on platforms without singularity
    #     "../envs/report_html.yaml"
    log:
        logs/"report/{capture_id}.metrics.log"
    script:
        "../scripts/report_capture_metrics.R"


rule report_capture:
    """
    Per-capture demux summary: cell counts, sample assignments, ambient RNA and
    agreement between the demultiplexing methods.
    """
    input:
        metrics=rules.report_capture_metrics.output.metrics,
        assignments=out_dir/"assignment/{capture_id}/sample_barcode_assignments.tsv",
        all_cells=out_dir/"assignment/{capture_id}/all_cells_stats.tsv",
        all_singlets=out_dir/"assignment/{capture_id}/all_singlets_stats.tsv",
        ambient_stats=out_dir/"vireo/{capture_id}/ambient_stats.csv",
        gt_distance=out_dir/"vireo/{capture_id}/fig_GT_distance_estimated.png",
        freemuxlet_corr=out_dir/"freemuxlet/{capture_id}/ref_clust_pearson_correlation.png",
        vireosnp_heatmap=out_dir/"souporcell/annotation/vireosnp/{capture_id}/matched_GPb_diff_scores_cluster_genotype.heatmap.png"
    output:
        # Single output: Snakemake permits an .Rmd script only one.
        report_dir/"captures/{capture_id}.html"
    params:
        pool_id=get_pool_for_capture
    container:
        config["containers"]["r"]
    # conda:  # swap for the container above on platforms without singularity
    #     "../envs/report_html.yaml"
    log:
        logs/"report/{capture_id}.log"
    script:
        "../scripts/report_capture.Rmd"


rule report_index:
    """
    Landing page: cross-capture overview plus tables searchable by capture and
    by sample.
    """
    input:
        metrics=expand(
            rules.report_capture_metrics.output.metrics,
            capture_id=captures.index.tolist()
        ),
        # The index links to every capture page, so depend on them: otherwise
        # only the index is reachable from `rule all` and the links 404.
        pages=expand(
            report_dir/"captures/{capture_id}.html",
            capture_id=captures.index.tolist()
        ),
        all_cells=out_dir/"assignment/report/all_cells_stats_full.tsv"
    output:
        report_dir/"index.html"
    container:
        config["containers"]["r"]
    # conda:  # swap for the container above on platforms without singularity
    #     "../envs/report_html.yaml"
    log:
        logs/"report/index.log"
    script:
        "../scripts/report_index.Rmd"
