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
    Return the pool_id wildcards.capture_id belongs to.
    """
    return captures.loc[wildcards.capture_id, "pool_id"]


rule report_capture:
    """
    Per-capture demux summary: cell counts, sample assignments, ambient RNA and
    agreement between the demultiplexing methods.
    """
    input:
        assignments=out_dir/"assignment/{capture_id}/sample_barcode_assignments.tsv",
        cell_assignment=out_dir/"assignment/{capture_id}/cell_assignment.tsv",
        all_cells=out_dir/"assignment/{capture_id}/all_cells_stats.tsv",
        all_singlets=out_dir/"assignment/{capture_id}/all_singlets_stats.tsv",
        ambient_stats=out_dir/"vireo/{capture_id}/ambient_stats.csv",
        ambient_summary=out_dir/"vireo/{capture_id}/ambient_summary.csv",
        souporcell_ambient=out_dir/"souporcell/{capture_id}/ambient_rna.txt",
        gt_distance=out_dir/"vireo/{capture_id}/fig_GT_distance_estimated.png",
        freemuxlet_corr=out_dir/"freemuxlet/{capture_id}/ref_clust_pearson_correlation.png",
        vireosnp_heatmap=out_dir/"souporcell/annotation/vireosnp/{capture_id}/matched_GPb_diff_scores_cluster_genotype.heatmap.png"
    output:
        html=report_dir/"captures/{capture_id}.html",
        # Small per-capture roll-up the index aggregates, so the index does not
        # have to re-read every capture's barcode-level table.
        metrics=report_dir/"metrics/{capture_id}.tsv"
    params:
        declared_samples=get_samples_in_capture,
        pool_id=get_pool_for_capture
    conda:
        "../envs/report_html.yaml"
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
            report_dir/"metrics/{capture_id}.tsv",
            capture_id=captures.index.tolist()
        ),
        all_cells=out_dir/"assignment/report/all_cells_stats_full.tsv"
    output:
        html=report_dir/"index.html"
    conda:
        "../envs/report_html.yaml"
    log:
        logs/"report/index.log"
    script:
        "../scripts/report_index.Rmd"
