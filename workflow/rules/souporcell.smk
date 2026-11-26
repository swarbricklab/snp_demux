# souporcell remap + vartrix (-> ref.mtx/alt.mtx/barcodes.tsv/variants.vcf) is
# produced upstream by sc_pileup (engine `souporcell`); this stage resumes at
# clustering. Pinned to cumulusprod/souporcell:2.5 to match the version snp_demux
# clustered with (keep the split behaviour-preserving; upgrade separately).
#
# --known_genotypes MUST be sc_pileup's variants.vcf (the covered subset the
# matrices were built on) so its rows align with alt.mtx/ref.mtx;
# --known_genotypes_sample_names restricts to the pool's donors.

SOUPORCELL_CONTAINER = config["containers"]["souporcell"]


rule souporcell_cluster:
    input:
        alt=pileup_path("souporcell", "alt.mtx"),
        ref=pileup_path("souporcell", "ref.mtx"),
        barcodes=pileup_path("souporcell", "barcodes.tsv"),
        variants=pileup_path("souporcell", "variants.vcf"),
    output:
        clusters_tmp=temp(out_dir / "souporcell/{capture_id}/clusters_tmp.tsv"),
    log:
        logs / "souporcell/{capture_id}/cluster.log",
    container:
        SOUPORCELL_CONTAINER
    params:
        k=get_number_of_samples_in_capture_pool,
        donors=get_donors_in_capture,
        min_ref=config["params"]["souporcell"]["min_ref"],
        min_alt=config["params"]["souporcell"]["min_alt"],
    shell:
        """
        # Restrict --known_genotypes_sample_names to donors actually present as
        # samples in variants.vcf. A donor absent from every genotype source has
        # no column here, and souporcell panics ("no entry found for key") on a
        # missing name. k is left unchanged: clusters beyond the supplied names
        # start at 0.5 at every site (no random start), so an ungenotyped donor
        # still gets a cluster.
        # NB: strict-mode-safe (snakemake wraps shell in `set -euo pipefail`).
        # `|| true` swallows the SIGPIPE from grep -m1 closing the reader; the
        # if/then keeps a not-found grep from aborting the loop.
        avail=$( (zcat -f {input.variants} 2>/dev/null || cat {input.variants}) | grep -m1 '^#CHROM' | cut -f10- | tr '\t' '\n') || true
        names=""
        for d in $(echo {params.donors} | tr ',' ' '); do
            if printf '%s\n' "$avail" | grep -qxF -- "$d"; then names="$names $d"; fi
        done
        # souporcell pairs --known_genotypes records with matrix rows by ORDER
        # (record i -> row i) and says nothing when the counts differ. sc_pileup
        # writes variants.vcf as the matrices' row list, so check it. Matrix
        # Market: the first line not starting with '%' is "rows cols entries".
        n_vcf=$( (zcat -f {input.variants} 2>/dev/null || cat {input.variants}) | grep -c -v -e '^#' -e '^$') || true
        n_alt=$(awk '!/^%/ {{ print $1; exit }}' {input.alt})
        n_ref=$(awk '!/^%/ {{ print $1; exit }}' {input.ref})
        if [ "$n_vcf" != "$n_alt" ] || [ "$n_ref" != "$n_alt" ]; then
            echo "ERROR: souporcell_cluster {wildcards.capture_id}: variants.vcf has $n_vcf records, alt.mtx $n_alt rows, ref.mtx $n_ref rows; --known_genotypes needs one record per matrix row" | tee {log} >&2
            exit 1
        fi
        echo "row check: $n_vcf variants.vcf records = $n_alt matrix rows" >{log}
        # One solve on one thread. With --known_genotypes souporcell starts each
        # cluster from a donor's genotypes (a cluster with no donor genotype at
        # 0.5) and never draws a random number
        # (init_cluster_centers_known_genotypes), so every restart is the same
        # solve and the first is kept (strict >). 100 restarts on 16 threads ran
        # 112 identical solves. --threads is a literal 1, not resources.threads:
        # each thread runs ceil(restarts/threads) solves. A run WITHOUT
        # --known_genotypes needs restarts, --seed and a fixed --threads again
        # (souporcell_cluster_known is one).
        /opt/souporcell/souporcell/target/release/souporcell -k {params.k} -a {input.alt} -r {input.ref} -b {input.barcodes} \
            --restarts 1 --min_ref {params.min_ref} --min_alt {params.min_alt} --threads 1 \
            --known_genotypes {input.variants} --known_genotypes_sample_names $names \
            >{output.clusters_tmp} 2>>{log}
        """


rule souporcell_troublet:
    input:
        alt=pileup_path("souporcell", "alt.mtx"),
        ref=pileup_path("souporcell", "ref.mtx"),
        clusters_tmp=rules.souporcell_cluster.output.clusters_tmp,
    output:
        clusters=out_dir / "souporcell/{capture_id}/clusters.tsv",
    log:
        logs / "souporcell/{capture_id}/troublet.log",
    container:
        SOUPORCELL_CONTAINER
    shell:
        """
        /opt/souporcell/troublet/target/release/troublet --alts {input.alt} --refs {input.ref} --clusters {input.clusters_tmp} \
            >{output.clusters} 2>{log}
        """


rule souporcell_consensus:
    input:
        alt=pileup_path("souporcell", "alt.mtx"),
        ref=pileup_path("souporcell", "ref.mtx"),
        variants=pileup_path("souporcell", "variants.vcf"),
        clusters=rules.souporcell_troublet.output.clusters,
    output:
        cluster_genotypes=out_dir / "souporcell/{capture_id}/cluster_genotypes.vcf",
        ambient_rna=out_dir / "souporcell/{capture_id}/ambient_rna.txt",
    log:
        logs / "souporcell/{capture_id}/consensus.log",
    container:
        SOUPORCELL_CONTAINER
    resources:
        # consensus.py: PBS memory used at most 7.1 GB over 17 captures (RSS
        # 2.0-6.0 GB on four of them); earlier union panels peaked ~11 GB and
        # were killed at 8 GB. Start at 8 GB and double per retry (retries: 2):
        # 8, 16, 32 GB; attempt 2 is the old flat 16 GB. A mem key for this rule
        # in profiles/workflow set-resources would override this on every
        # attempt: keep it out.
        mem_mb=lambda wildcards, attempt: 8000 * 2 ** (attempt - 1),
    params:
        ploidy=config["params"]["souporcell"]["ploidy"],
    shell:
        """
        outdir=$(dirname {output.cluster_genotypes})
        /opt/souporcell/consensus.py -c {input.clusters} -a {input.alt} -r {input.ref} -p {params.ploidy} \
            --output_dir $outdir --soup_out {output.ambient_rna} \
            --vcf_out {output.cluster_genotypes} --vcf {input.variants} \
            >{log} 2>&1
        """


# ---------------------------------------------------------------------------
# souporcell second run, DE NOVO (parallel to the rules above).
# ---------------------------------------------------------------------------
# The rules above start from genotypes: sc_pileup's variants.vcf, passed to
# --known_genotypes, carries the donor columns. This run passes no genotypes,
# so it clusters de novo from random starts, and hierarchical_consensus
# compares the two partitions (denovo_known_agree, the only output it feeds).
# Outputs keep the souporcell/{c}/known/ path, and the column keeps its name.
rule souporcell_cluster_known:
    input:
        alt=pileup_path("souporcell", "alt.mtx"),
        ref=pileup_path("souporcell", "ref.mtx"),
        barcodes=pileup_path("souporcell", "barcodes.tsv"),
    output:
        clusters_tmp=temp(out_dir / "souporcell/{capture_id}/known/clusters_tmp.tsv"),
    log:
        logs / "souporcell/{capture_id}/known/cluster.log",
    container:
        SOUPORCELL_CONTAINER
    resources:
        threads=SOUPORCELL_DENOVO_THREADS,
    params:
        k=get_number_of_samples_in_capture_pool,
        min_ref=config["params"]["souporcell"]["min_ref"],
        min_alt=config["params"]["souporcell"]["min_alt"],
        restarts=SOUPORCELL_DENOVO_RESTARTS,
        seed=SOUPORCELL_SEED,
        threads=SOUPORCELL_DENOVO_THREADS,
    shell:
        """
        # De novo: the result depends on --restarts, --seed and --threads (each
        # thread draws its own starts from the seed and runs
        # ceil(restarts/threads) solves). All three come from config
        # (params.souporcell restarts, seed, denovo_threads), not from the
        # profile; the CPU request follows denovo_threads through resources.
        /opt/souporcell/souporcell/target/release/souporcell -k {params.k} -a {input.alt} -r {input.ref} -b {input.barcodes} \
            --restarts {params.restarts} --min_ref {params.min_ref} --min_alt {params.min_alt} \
            --threads {params.threads} --seed {params.seed} \
            >{output.clusters_tmp} 2>{log}
        """


rule souporcell_troublet_known:
    input:
        alt=pileup_path("souporcell", "alt.mtx"),
        ref=pileup_path("souporcell", "ref.mtx"),
        clusters_tmp=rules.souporcell_cluster_known.output.clusters_tmp,
    output:
        clusters=out_dir / "souporcell/{capture_id}/known/clusters.tsv",
    log:
        logs / "souporcell/{capture_id}/known/troublet.log",
    container:
        SOUPORCELL_CONTAINER
    shell:
        """
        /opt/souporcell/troublet/target/release/troublet --alts {input.alt} --refs {input.ref} --clusters {input.clusters_tmp} \
            >{output.clusters} 2>{log}
        """


rule souporcell_consensus_known:
    input:
        alt=pileup_path("souporcell", "alt.mtx"),
        ref=pileup_path("souporcell", "ref.mtx"),
        variants=pileup_path("souporcell", "variants.vcf"),
        clusters=rules.souporcell_troublet_known.output.clusters,
    output:
        cluster_genotypes=out_dir
        / "souporcell/{capture_id}/known/cluster_genotypes.vcf",
        ambient_rna=out_dir / "souporcell/{capture_id}/known/ambient_rna.txt",
    log:
        logs / "souporcell/{capture_id}/known/consensus.log",
    container:
        SOUPORCELL_CONTAINER
    resources:
        # As souporcell_consensus: 8, 16, 32 GB by attempt; keep mem out of the
        # profile. Nothing in rule all requests this rule's outputs; it is kept
        # consistent with souporcell_consensus.
        mem_mb=lambda wildcards, attempt: 8000 * 2 ** (attempt - 1),
    params:
        ploidy=config["params"]["souporcell"]["ploidy"],
    shell:
        """
        outdir=$(dirname {output.cluster_genotypes})
        /opt/souporcell/consensus.py -c {input.clusters} -a {input.alt} -r {input.ref} -p {params.ploidy} \
            --output_dir $outdir --soup_out {output.ambient_rna} \
            --vcf_out {output.cluster_genotypes} --vcf {input.variants} \
            >{log} 2>&1
        """
