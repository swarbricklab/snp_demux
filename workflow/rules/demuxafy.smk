# ==================================
# demuxafy
# ----------------------------------
# Run demuxafy "combine results" to find common cell annotation (to individuals) across methods
# ==================================
rule demuxafy:
    input:
        demuxlet=rules.demuxlet.output.best,
        freemuxlet=rules.freemuxlet.output.samples,
        souporcell=rules.souporcell.output.clusters,
        vireo=rules.vireo.output.summary
    output:
        # TODO: add wildcard based logic to determine the combine method (e.g., AtLeastHalfSinglet)
        combined_results=temp(out_dir/"assignment/{capture_id}/combined_results.tsv"),
        combined_resultsSinglets_upset=temp(out_dir/"assignment/{capture_id}/combined_resultsSinglets_upset.pdf"),
        combined_results_annotated=temp(out_dir/"assignment/{capture_id}/combined_results_w_combined_assignments.tsv"),
        combined_results_demultiplexing_summary=temp(out_dir/"assignment/{capture_id}/combined_results_demultiplexing_summary.tsv"),
        combined_results_summary=temp(out_dir/"assignment/{capture_id}/combined_results_summary.tsv")
    params:
        freemuxlet_correlation_limit=config['params']['demuxafy']['freemuxlet_correlation_limit'],
        souporcell_correlation_limit=config['params']['demuxafy']['souporcell_correlation_limit'],
        pct_agreement=config['params']['demuxafy']['pct_agreement']
    container:
        config["containers"]["demuxafy"]
    log:
        logs/"assignment/{capture_id}/demuxafy.log"
    shell:
        """
        demuxlet_dir=$(dirname {input.demuxlet})
        freemuxlet_dir=$(dirname {input.freemuxlet})
        souporcell_dir=$(dirname {input.souporcell})
        vireo_dir=$(dirname {input.vireo})
        Combine_Results.R \
            -o {output.combined_results} \
            --demuxlet $demuxlet_dir \
            --freemuxlet $freemuxlet_dir \
            --freemuxlet_correlation_limit {params.freemuxlet_correlation_limit} \
            --souporcell $souporcell_dir \
            --souporcell_correlation_limit {params.souporcell_correlation_limit} \
            --vireo $vireo_dir \
            --method "AtLeastHalfSinglet" \
            --pct_agreement {params.pct_agreement} \
            > {log}
        """

# Convert PDF plot above to PNG format
rule pdf_to_png:
    input:
        Singlets_upset_pdf=rules.demuxafy.output.combined_resultsSinglets_upset
    output:
        Singlets_upset_png=out_dir/"assignment/{capture_id}/plots/combined_resultsSinglets_upset.png"
    conda:
        "../envs/pdftools.yaml"
    log:
        logs/"pdf_to_png/{capture_id}.log"
    script:
        "../scripts/pdf_to_png.R"


# ==================================
# add_vireosnp_to_combined_results
# ----------------------------------
# Add VireoSNP annotation results to the demuxafy combined annotations
# ==================================
rule add_vireosnp_to_combined_results:
    input:
        demuxafy_assignments=rules.demuxafy.output.combined_results_annotated,
        vireosnp_assignments=rules.vireosnp.output.clusters_annotated
    output:
        assignments_with_vireosnp=temp(out_dir/"assignment/{capture_id}/combined_results_w_combined_assignments.inc_vireosnp.tsv")
    conda:
        "../envs/report.yaml"
    log:
        logs/"assignment/{capture_id}/add_vireosnp_to_combined_results.log"
    script:
        "../scripts/add_vireosnp_to_combined_results.R"


# # ==================================
# # wrangle_demuxafy_output
# # ----------------------------------
# # post-processing of demuxafy results and summary stats from cells and singlets
# # ==================================
rule wrangle_demuxafy_output:
    input:
        metadata=config["deps"]["samples"],
        demuxafy_assignments=rules.add_vireosnp_to_combined_results.output.assignments_with_vireosnp
    output:
        # TODO: think about more readable/informative name
        combined_results_w_combined_assignments_edited=out_dir/"assignment/{capture_id}/sample_barcode_assignments.tsv",
        all_cells_stats=out_dir/"assignment/{capture_id}/all_cells_stats.tsv",
        all_singlets_stats=out_dir/"assignment/{capture_id}/all_singlets_stats.tsv",
        cell_assignment=out_dir/"assignment/{capture_id}/cell_assignment.tsv"
    params:
        # TODO: add wildcard based processing to demuxafy "method"
        assignment_column="AtLeastHalfSinglet_Individual_Assignment",
        droplettype_column="AtLeastHalfSinglet_DropletType",
        spikein_sample_id=config["spikein_sample_id"]
    conda:
        "../envs/report.yaml"
    log:
        logs/"assignment/{capture_id}/wrangle_demuxafy_output.log"
    script:
        "../scripts/wrangle_demuxafy_output.R"

rule demuxafy_count_cells_AtLeastHalfSinglet:
    input:
        clusters_annotated=rules.wrangle_demuxafy_output.output.combined_results_w_combined_assignments_edited
    output:
        cell_counts_all=out_dir/"assignment/{capture_id}/plots/demuxafy_cell_count_all.png",
        cell_counts_singlets=out_dir/"assignment/{capture_id}/plots/demuxafy_cell_count_singlets.png",
        table_cell_counts_all=out_dir/"assignment/{capture_id}/demuxafy_cell_count_table.tsv"
    params:
        dir_output=out_dir/"assignment/{capture_id}/"
    conda:
        "../envs/vireosnp.yaml"
    log:
        logs/"assignment/{capture_id}/output.count_cells.log"
    script:
        "../scripts/demuxafy_counts_AtLeastHalfSinglet.py"