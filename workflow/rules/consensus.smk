# ==================================
# consensus
# ----------------------------------
# Per-barcode consensus of demuxlet, freemuxlet, souporcell and vireo (plus the
# vireoSNP annotation of souporcell's clusters, which does not vote). Replaces
# Demuxafy's Combine_Results.R; the rule is params.consensus.method (see
# workflow/scripts/consensus.py).
# ==================================
CONSENSUS = {
    "method": "informative_majority",
    "min_informative": 2,
    "cluster_naming": "genotype",
    "min_genotype_margin": 0.1,
    "recover_from": ["souporcell", "vireo"],
    **config["params"].get("consensus", {}),
}


rule consensus:
    input:
        demuxlet=rules.demuxlet.output.best,
        freemuxlet=rules.freemuxlet.output.samples,
        souporcell=rules.souporcell_troublet.output.clusters,
        vireo=rules.relabel_vireo.output.donor_ids,
        vireosnp=rules.vireosnp.output.clusters_annotated,
        souporcell_scores=rules.vireosnp.output.scores_full,
        freemuxlet_scores=rules.freemuxlet_genotype_match.output.scores,
        genotyped_donors=rules.extract_genotyped_donors.output.genotyped_donors,
    output:
        calls=temp(out_dir / "assignment/{capture_id}/consensus_calls.tsv"),
    log:
        logs / "assignment/{capture_id}/consensus.log",
    container:
        config["containers"]["vireosnp"]
    params:
        pct_agreement=config["params"]["demuxafy"]["pct_agreement"],
        method=CONSENSUS["method"],
        min_informative=CONSENSUS["min_informative"],
        cluster_naming=CONSENSUS["cluster_naming"],
        min_genotype_margin=CONSENSUS["min_genotype_margin"],
        recover_from=CONSENSUS["recover_from"],
        pool_donors=get_donors_in_capture,
    script:
        "../scripts/consensus.py"


# Append vireo's and souporcell's per-barcode posteriors (in donor space, so
# wrangle_assignments.R recodes their best-singlet columns to sample ids with
# the rest). See workflow/scripts/add_posteriors.py.
rule add_posteriors:
    input:
        table=rules.consensus.output.calls,
        vireo=rules.relabel_vireo.output.donor_ids,
        souporcell=rules.souporcell_troublet.output.clusters,
        samples=config["deps"]["samples"],
    output:
        table=temp(
            out_dir / "assignment/{capture_id}/consensus_calls.inc_posteriors.tsv"
        ),
    log:
        logs / "assignment/{capture_id}/add_posteriors.log",
    container:
        config["containers"]["vireosnp"]
    script:
        "../scripts/add_posteriors.py"


rule wrangle_demuxafy_assignments:
    input:
        sample_map=config["deps"]["samples"],
        #vcf=config["deps"]["vcf"],
        demuxafy_assignments=rules.add_posteriors.output.table,
        # With params.elimination_gate.action: withhold, a FAILed by-elimination
        # donor's cells are reported as unassigned in cell_assignment.tsv.
        elimination_gate=rules.elimination_gate.output.tsv,
    output:
        sample_barcode_assignments=out_dir
        / "assignment/{capture_id}/sample_barcode_assignments.tsv",
        cell_assignment=out_dir / "assignment/{capture_id}/cell_assignment.tsv",
    # conda:  # conda fallback for platforms without singularity (use instead of the container)
    #     "../envs/wrangle_demuxafy.yaml"
    log:
        logs / "assignment/{capture_id}/wrangle_demuxafy_output.assignments.log",
    container:
        config["containers"]["r"]
    params:
        assignment_column="Consensus_Individual_Assignment",
        droplettype_column="Consensus_DropletType",
    script:
        "../scripts/wrangle_assignments.R"


rule wrangle_stats:
    input:
        sample_barcode_assignments=rules.wrangle_demuxafy_assignments.output.sample_barcode_assignments,
    output:
        all_cells_stats=out_dir / "assignment/{capture_id}/all_cells_stats.tsv",
        all_singlets_stats=out_dir / "assignment/{capture_id}/all_singlets_stats.tsv",
    # conda:  # conda fallback for platforms without singularity (use instead of the container)
    #     "../envs/wrangle_demuxafy.yaml"
    log:
        logs / "assignment/{capture_id}/wrangle_stats.log",
    container:
        config["containers"]["r"]
    params:
        assignment_column="Consensus_Individual_Assignment",
        droplettype_column="Consensus_DropletType",
    script:
        "../scripts/wrangle_stats.R"


rule demuxafy_count_cells_AtLeastHalfSinglet:
    input:
        clusters_annotated=rules.wrangle_demuxafy_assignments.output.sample_barcode_assignments,
    output:
        cell_counts_all=out_dir
        / "assignment/{capture_id}/plots/demuxafy_cell_count_all.png",
        cell_counts_singlets=out_dir
        / "assignment/{capture_id}/plots/demuxafy_cell_count_singlets.png",
        table_cell_counts_all=out_dir
        / "assignment/{capture_id}/demuxafy_cell_count_table.tsv",
    # conda:  # conda fallback for platforms without singularity (use instead of the container)
    #     "../envs/vireosnp.yaml"
    log:
        logs / "assignment/{capture_id}/output.count_cells.log",
    container:
        config["containers"]["vireosnp"]
    params:
        dir_output=out_dir / "assignment/{capture_id}/",
    script:
        "../scripts/demuxafy_counts_AtLeastHalfSinglet.py"
