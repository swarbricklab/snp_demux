# ==================================
# hierarchical_consensus
# ----------------------------------
# Sex-aware hierarchical consensus, run ALONGSIDE demuxafy (demuxafy.smk is left
# untouched so the PR can show both side by side). De-novo souporcell decides
# cluster membership (recovery); vireosnp GPb_diff margins decide the donor label;
# panel-derived donor sex breaks low-margin ties. See
# docs/hierarchical_demux_design.md.
# ==================================

# Join the per-capture sex_calls.csv onto the de-novo souporcell barcodes. The sex
# CSV path comes via params (sex_path()) not input: it is an optional upstream
# source with no producing rule, and join_sex.py degrades to empty sex columns if
# it is absent.
rule join_sex:
    group: "postproc"
    input:
        base_table=rules.souporcell_troublet.output.clusters,
    output:
        sex_joined=out_dir/"assignment/{capture_id}/sex_calls.joined.tsv",
    params:
        sex_csv=sex_path(),
    conda:
        "../envs/vireosnp.yaml"
    log:
        logs/"assignment/{capture_id}/join_sex.log"
    script:
        "../scripts/join_sex.py"


rule hierarchical_consensus:
    group: "postproc"
    input:
        # recovery backbone (de-novo souporcell) + doublet/singlet log-probs
        denovo_clusters=rules.souporcell_troublet.output.clusters,
        # known-GT souporcell run -> de-novo-vs-known partition agreement QC
        known_clusters=rules.souporcell_troublet_known.output.clusters,
        # vireosnp cluster<->donor annotation + GPb_diff margin matrix
        vireosnp_annotated=rules.vireosnp.output.clusters_annotated,
        scores_full=rules.vireosnp.output.scores_full,
        # vireo per-cell doublet posterior
        vireo_donor_ids=rules.relabel_vireo.output.donor_ids,
        # real donor genotypes -> donor sex (chrX het / chrY presence)
        known_genotypes_vcf=rules.known_genotypes.output.vcf,
        # per-barcode sex join for this capture
        sex_joined=rules.join_sex.output.sex_joined,
    output:
        hier=out_dir/"assignment/{capture_id}/hierarchical_consensus.tsv",
    params:
        margin_min=config['params']['hierarchical']['margin_min'],
        sex_weight=config['params']['hierarchical']['sex_weight'],
        min_cells_floor=config['params']['hierarchical']['min_cells_floor'],
        low_conf_keep=config['params']['hierarchical']['low_conf_keep'],
    conda:
        "../envs/vireosnp.yaml"
    log:
        logs/"assignment/{capture_id}/hierarchical_consensus.log"
    script:
        "../scripts/hierarchical_consensus.py"
