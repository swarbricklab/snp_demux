# ==================================
# Beata Kiedik
# ----------------------------------
# CREATED: 14/12/2022
# ----------------------------------
# Rules for counting cells in SNP clusters, based on known genotyping
# ------------------------------------------------------------------------------
# ==================================

rule demuxafy_count_cells_AtLeastHalfSinglet:
    input:
        clusters_annotated=out_dir/"assignment/{capture_id}/combined_results_w_combined_assignments_edited.tsv"
    output:
        cell_counts_all=out_dir/"assignment/{capture_id}/demuxafy_cell_count_all.png",
        cell_counts_singlets=out_dir/"assignment/{capture_id}/demuxafy_cell_count_singlets.png",
        table_cell_counts_all=out_dir/"assignment/{capture_id}/demuxafy_cell_count_table.tsv"
    params:
        dir_output=out_dir/"assignment/{capture_id}/"
    conda:
        "../envs/vireosnp.yaml"
    log:
        logs/"assignment/{capture_id}/output.count_cells.log"
    script:
        "../scripts/demuxafy_counts_AtLeastHalfSinglet.py"
