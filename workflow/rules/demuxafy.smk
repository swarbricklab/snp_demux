# ==================================
# demuxafy
# ----------------------------------
# Run demuxafy "combine results" to find common cell annotation (to individuals) across methods
# ==================================
# Per-tool outputs handed to Demuxafy's Combine_Results.R. Kept as a single
# source of truth so the input files (below) and the shell --<tool> <dir> flags
# stay in lock-step. NB: Combine_Results.R takes a DIRECTORY per tool (the rule
# passes dirname(<file>)); flag names are literally --demuxlet/--freemuxlet/
# --souporcell/--vireo.
_DEMUXAFY_TOOL_OUTPUTS = {
    "demuxlet": rules.demuxlet.output.best,
    "freemuxlet": rules.freemuxlet.output.samples,
    "souporcell": rules.souporcell_troublet.output.clusters,
    "vireo": rules.relabel_vireo.output.donor_ids,
}


def _demuxafy_require_two_tools():
    """Demuxafy consensus needs >=2 per-cell tools to combine. With fewer, it
    cannot form a consensus; fail early with a clear message rather than crashing
    inside Combine_Results.R. (Default is 4 tools, so this never triggers in the
    backward-compatible case.)"""
    if len(enabled_tools()) < 2:
        raise WorkflowError(
            "Demuxafy consensus (rule demuxafy) requires at least 2 enabled "
            "per-cell demux tools in config params.tools, but only %s "
            "is/are enabled. Enable a second tool or drop the demuxafy targets."
            % enabled_tools()
        )


def demuxafy_inputs(wildcards):
    """Only the ENABLED tools' outputs are declared as inputs, so disabled tools
    (e.g. demuxlet/freemuxlet without popscle pileups) are never demanded."""
    _demuxafy_require_two_tools()
    c = wildcards.capture_id
    return {
        tool: _DEMUXAFY_TOOL_OUTPUTS[tool].format(capture_id=c)
        for tool in enabled_tools()
    }


def demuxafy_combine_flags(wildcards):
    """Build the per-tool Combine_Results.R flags for just the enabled tools,
    preserving the original flag order/values so the invocation is byte-for-byte
    identical when all four tools are enabled."""
    _demuxafy_require_two_tools()
    c = wildcards.capture_id
    p = config['params']['demuxafy']
    parts = []
    if tool_enabled("demuxlet"):
        parts.append("--demuxlet " + str(Path(_DEMUXAFY_TOOL_OUTPUTS["demuxlet"].format(capture_id=c)).parent))
    if tool_enabled("freemuxlet"):
        parts.append("--freemuxlet " + str(Path(_DEMUXAFY_TOOL_OUTPUTS["freemuxlet"].format(capture_id=c)).parent))
        parts.append("--freemuxlet_correlation_limit " + str(p['freemuxlet_correlation_limit']))
    if tool_enabled("souporcell"):
        parts.append("--souporcell " + str(Path(_DEMUXAFY_TOOL_OUTPUTS["souporcell"].format(capture_id=c)).parent))
        parts.append("--souporcell_correlation_limit " + str(p['souporcell_correlation_limit']))
    if tool_enabled("vireo"):
        parts.append("--vireo " + str(Path(_DEMUXAFY_TOOL_OUTPUTS["vireo"].format(capture_id=c)).parent))
    return " ".join(parts)


def demuxafy_ref():
    """Combine_Results.R --ref must name an ENABLED tool. Prefer Vireo, else
    Souporcell, else the first enabled tool (capitalised, e.g. Demuxlet)."""
    if tool_enabled("vireo"):
        return "Vireo"
    if tool_enabled("souporcell"):
        return "Souporcell"
    return enabled_tools()[0].capitalize()


rule demuxafy:
    group: "postproc"
    input:
        unpack(demuxafy_inputs)
    output:
        # TODO: add wildcard based logic to determine the combine method (e.g., AtLeastHalfSinglet)
        combined_results=temp(out_dir/"assignment/{capture_id}/combined_results.tsv"),
        combined_resultsSinglets_upset=temp(out_dir/"assignment/{capture_id}/combined_resultsSinglets_upset.pdf"),
        combined_results_annotated=temp(out_dir/"assignment/{capture_id}/combined_results_w_combined_assignments.tsv"),
        combined_results_demultiplexing_summary=temp(out_dir/"assignment/{capture_id}/combined_results_demultiplexing_summary.tsv"),
        combined_results_summary=temp(out_dir/"assignment/{capture_id}/combined_results_summary.tsv")
    params:
        combine_flags=demuxafy_combine_flags,
        ref=demuxafy_ref(),
        pct_agreement=config['params']['demuxafy']['pct_agreement']
    container:
        config["containers"]["demuxafy"]
    log:
        logs/"assignment/{capture_id}/demuxafy.log"
    shell:
        """
        Combine_Results.R \
            -o {output.combined_results} \
            {params.combine_flags} \
            --method "AtLeastHalfSinglet" \
            --ref {params.ref} \
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
    group: "postproc"
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

rule wrangle_demuxafy_assignments:
    group: "postproc"
    input:
        sample_map=config["deps"]["samples"],
        #vcf=config["deps"]["vcf"],
        demuxafy_assignments=rules.add_vireosnp_to_combined_results.output.assignments_with_vireosnp
    output:
        sample_barcode_assignments=out_dir/"assignment/{capture_id}/sample_barcode_assignments.tsv",
        cell_assignment=out_dir/"assignment/{capture_id}/cell_assignment.tsv"
    params:
        assignment_column="AtLeastHalfSinglet_Individual_Assignment",
        droplettype_column="AtLeastHalfSinglet_DropletType",
    conda:
        "../envs/wrangle_demuxafy.yaml"
    log:
        logs/"assignment/{capture_id}/wrangle_demuxafy_output.assignments.log"
    script:
        "../scripts/wrangle_assignments.R"

rule wrangle_stats:
    group: "postproc"
    input:
        sample_barcode_assignments=rules.wrangle_demuxafy_assignments.output.sample_barcode_assignments
    output:
        all_cells_stats=out_dir/"assignment/{capture_id}/all_cells_stats.tsv",
        all_singlets_stats=out_dir/"assignment/{capture_id}/all_singlets_stats.tsv"
    params:
        assignment_column="AtLeastHalfSinglet_Individual_Assignment"
    conda:
        "../envs/wrangle_demuxafy.yaml"
    log:
        logs/"assignment/{capture_id}/wrangle_stats.log"
    script:
        "../scripts/wrangle_stats.R"

rule demuxafy_count_cells_AtLeastHalfSinglet:
    input:
        clusters_annotated=rules.wrangle_demuxafy_assignments.output.sample_barcode_assignments
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