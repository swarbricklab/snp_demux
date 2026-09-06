# souporcell remap + vartrix (-> ref.mtx/alt.mtx/barcodes.tsv/variants.vcf) is
# produced upstream by sc_pileup (engine `souporcell`); this stage resumes at
# clustering. Pinned to cumulusprod/souporcell:2.5 to match the version snp_demux
# clustered with (keep the split behaviour-preserving; upgrade separately).
#
# --known_genotypes MUST be sc_pileup's variants.vcf (the covered subset the
# matrices were built on) so its rows align with alt.mtx/ref.mtx;
# --known_genotypes_sample_names restricts to the pool's donors.

SOUPORCELL_CONTAINER = "docker://cumulusprod/souporcell:2.5"


rule souporcell_cluster:
    input:
        alt=pileup_path("souporcell", "alt.mtx"),
        ref=pileup_path("souporcell", "ref.mtx"),
        barcodes=pileup_path("souporcell", "barcodes.tsv"),
        variants=pileup_path("souporcell", "variants.vcf"),
    output:
        clusters_tmp=temp(out_dir/"souporcell/{capture_id}/clusters_tmp.tsv"),
    params:
        k=get_number_of_samples_in_capture_pool,
        donors=get_donors_in_capture,
        restarts=config['params']['souporcell']['restarts'],
        min_ref=config['params']['souporcell']['min_ref'],
        min_alt=config['params']['souporcell']['min_alt'],
    container:
        SOUPORCELL_CONTAINER
    log:
        logs/"souporcell/{capture_id}/cluster.log"
    shell:
        """
        names=$(echo {params.donors} | tr ',' ' ')
        /opt/souporcell/souporcell/target/release/souporcell -k {params.k} -a {input.alt} -r {input.ref} -b {input.barcodes} \
            --restarts {params.restarts} --min_ref {params.min_ref} --min_alt {params.min_alt} \
            --known_genotypes {input.variants} --known_genotypes_sample_names $names \
            > {output.clusters_tmp} 2> {log}
        """


rule souporcell_troublet:
    input:
        alt=pileup_path("souporcell", "alt.mtx"),
        ref=pileup_path("souporcell", "ref.mtx"),
        clusters_tmp=rules.souporcell_cluster.output.clusters_tmp,
    output:
        clusters=out_dir/"souporcell/{capture_id}/clusters.tsv",
    container:
        SOUPORCELL_CONTAINER
    log:
        logs/"souporcell/{capture_id}/troublet.log"
    shell:
        """
        /opt/souporcell/troublet/target/release/troublet --alts {input.alt} --refs {input.ref} --clusters {input.clusters_tmp} \
            > {output.clusters} 2> {log}
        """


rule souporcell_consensus:
    input:
        alt=pileup_path("souporcell", "alt.mtx"),
        ref=pileup_path("souporcell", "ref.mtx"),
        variants=pileup_path("souporcell", "variants.vcf"),
        clusters=rules.souporcell_troublet.output.clusters,
    output:
        cluster_genotypes=out_dir/"souporcell/{capture_id}/cluster_genotypes.vcf",
        ambient_rna=out_dir/"souporcell/{capture_id}/ambient_rna.txt",
    params:
        ploidy=config['params']['souporcell']['ploidy'],
    container:
        SOUPORCELL_CONTAINER
    log:
        logs/"souporcell/{capture_id}/consensus.log"
    shell:
        """
        outdir=$(dirname {output.cluster_genotypes})
        /opt/souporcell/consensus.py -c {input.clusters} -a {input.alt} -r {input.ref} -p {params.ploidy} \
            --output_dir $outdir --soup_out {output.ambient_rna} \
            --vcf_out {output.cluster_genotypes} --vcf {input.variants} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# souporcell KNOWN-GT second mode (parallel to the de-novo rules above).
# ---------------------------------------------------------------------------
# The de-novo rules above pass sc_pileup's variants.vcf (a bare site list) to
# --known_genotypes, so clustering is genotype-free (the recovery backbone).
# These parallel rules instead point --known_genotypes at the REAL per-pool donor
# genotypes (rule known_genotypes in split.smk = the informative panel subset to
# the pool's donors), giving a direct genotype-anchored clustering to compare
# against de-novo. Outputs live under souporcell/{c}/known/ so the de-novo path
# is untouched. NB: known_genotypes.vcf is built from the SAME informative panel
# sc_pileup piled against (see panel_path docstring), so its sites match the
# alt/ref matrices by construction; souporcell matches genotype rows by position.
rule souporcell_cluster_known:
    input:
        alt=pileup_path("souporcell", "alt.mtx"),
        ref=pileup_path("souporcell", "ref.mtx"),
        barcodes=pileup_path("souporcell", "barcodes.tsv"),
        known_gt=rules.known_genotypes.output.vcf,
    output:
        clusters_tmp=temp(out_dir/"souporcell/{capture_id}/known/clusters_tmp.tsv"),
    params:
        k=get_number_of_samples_in_capture_pool,
        donors=get_donors_in_capture,
        restarts=config['params']['souporcell']['restarts'],
        min_ref=config['params']['souporcell']['min_ref'],
        min_alt=config['params']['souporcell']['min_alt'],
    container:
        SOUPORCELL_CONTAINER
    log:
        logs/"souporcell/{capture_id}/known/cluster.log"
    shell:
        """
        names=$(echo {params.donors} | tr ',' ' ')
        /opt/souporcell/souporcell/target/release/souporcell -k {params.k} -a {input.alt} -r {input.ref} -b {input.barcodes} \
            --restarts {params.restarts} --min_ref {params.min_ref} --min_alt {params.min_alt} \
            --known_genotypes {input.known_gt} --known_genotypes_sample_names $names \
            > {output.clusters_tmp} 2> {log}
        """


rule souporcell_troublet_known:
    input:
        alt=pileup_path("souporcell", "alt.mtx"),
        ref=pileup_path("souporcell", "ref.mtx"),
        clusters_tmp=rules.souporcell_cluster_known.output.clusters_tmp,
    output:
        clusters=out_dir/"souporcell/{capture_id}/known/clusters.tsv",
    container:
        SOUPORCELL_CONTAINER
    log:
        logs/"souporcell/{capture_id}/known/troublet.log"
    shell:
        """
        /opt/souporcell/troublet/target/release/troublet --alts {input.alt} --refs {input.ref} --clusters {input.clusters_tmp} \
            > {output.clusters} 2> {log}
        """


rule souporcell_consensus_known:
    input:
        alt=pileup_path("souporcell", "alt.mtx"),
        ref=pileup_path("souporcell", "ref.mtx"),
        variants=pileup_path("souporcell", "variants.vcf"),
        clusters=rules.souporcell_troublet_known.output.clusters,
    output:
        cluster_genotypes=out_dir/"souporcell/{capture_id}/known/cluster_genotypes.vcf",
        ambient_rna=out_dir/"souporcell/{capture_id}/known/ambient_rna.txt",
    params:
        ploidy=config['params']['souporcell']['ploidy'],
    container:
        SOUPORCELL_CONTAINER
    log:
        logs/"souporcell/{capture_id}/known/consensus.log"
    shell:
        """
        outdir=$(dirname {output.cluster_genotypes})
        /opt/souporcell/consensus.py -c {input.clusters} -a {input.alt} -r {input.ref} -p {params.ploidy} \
            --output_dir $outdir --soup_out {output.ambient_rna} \
            --vcf_out {output.cluster_genotypes} --vcf {input.variants} \
            > {log} 2>&1
        """
