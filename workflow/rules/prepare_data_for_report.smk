# ==================================
# Beata Kiedik 
# ----------------------------------
# CREATED: 07/03/2023
# ----------------------------------
# Rules for preparing and combining stats for summary report
# ------------------------------------------------------------------------------
# ==================================

# ==================================
# combine cell statistics across captures
# ----------------------------------
# NOTE/TODO: this rule logic and script has been re-used from the sc-workflow - think about best way to re-use/centralise these (maybe as a wrapper rule/script)
# ==================================
rule combine_cell_assignments:
    input:
        data_files=expand(out_dir/"assignment/{capture_id}/sample_barcode_assignments.tsv", capture_id=captures.index.tolist())
    output:
        # TODO: should this be defined as an input to the report rule/script?
        combined=out_dir/"assignment/report/combined_results.tsv"
    params:
        # Positionally paired with input.data_files -- expand() below iterates
        # this same list, in this order.
        capture_ids=captures.index.tolist()
    log:
        logs/"assignment/report/output.combine_cell_assignments.log"
    container:
        config["containers"]["r"]
    # conda:  # swap for the container above on platforms without singularity
    #     "../envs/combine_matrices.yaml"
    script:
        "../scripts/combine_matrices_by_row.R"

use rule combine_cell_assignments as combine_cell_stats with:
    input:
        data_files=expand(out_dir/"assignment/{capture_id}/all_cells_stats.tsv", capture_id=captures.index.tolist())
    output:
        combined=out_dir/"assignment/report/all_cells_stats_full.tsv"
    log:
        logs/"assignment/report/output.all_cells_stats_full.tsv.log"

use rule combine_cell_assignments as combine_singlet_stats with:
    input:
        data_files=expand(out_dir/"assignment/{capture_id}/all_singlets_stats.tsv", capture_id=captures.index.tolist())
    output:
        combined=out_dir/"assignment/report/all_singlets_stats_full.tsv"
    log:
        logs/"assignment/report/output.all_singlet_stats_full.tsv.log"
