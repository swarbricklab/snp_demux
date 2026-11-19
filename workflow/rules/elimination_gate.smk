# ==================================
# Concordance gate for the by-elimination missing-donor assignment
# ----------------------------------
# When one pool donor has no reference genotype, vireo learns it and the
# souporcell annotation gives it the leftover cluster, both by elimination.
# This gate checks that the leftover cluster's genotype actually agrees with the
# learned one before the assignment is trusted. See scripts/elimination_gate.py.
# ==================================


rule elimination_gate:
    input:
        genotyped_donors=rules.extract_genotyped_donors.output.genotyped_donors,
        mapping=rules.vireosnp.output.final_sample_cluster_mappings,
        vireo_vcf=rules.reheader_vireo_vcf.output.vcf,
        cluster_genotypes=rules.souporcell_consensus.output.cluster_genotypes,
    output:
        tsv=out_dir / "assignment/{capture_id}/elimination_gate.tsv",
    log:
        logs / "assignment/{capture_id}/elimination_gate.log",
    container:
        config["containers"]["vireosnp"]
    params:
        donors=get_donors_in_capture,
        min_concordance=elimination_gate_param("min_concordance"),
        max_drop=elimination_gate_param("max_drop"),
        min_sites=elimination_gate_param("min_sites"),
        action=elimination_gate_param("action"),
    script:
        "../scripts/elimination_gate.py"


rule combine_elimination_gate:
    input:
        gates=expand(
            out_dir / "assignment/{capture_id}/elimination_gate.tsv",
            capture_id=captures.index.tolist(),
        ),
        # Positionally paired with params.capture_ids, for cross-capture
        # stability of a missing donor's learned genotype.
        vireo=expand(
            out_dir / "vireo/{capture_id}/GT_donors.vireo.vcf.gz",
            capture_id=captures.index.tolist(),
        ),
    output:
        tsv=out_dir / "assignment/report/elimination_gate.tsv",
    log:
        logs / "assignment/report/elimination_gate.log",
    container:
        config["containers"]["vireosnp"]
    params:
        capture_ids=captures.index.tolist(),
    script:
        "../scripts/elimination_gate.py"
