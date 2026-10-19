# popscle dsc-pileup is produced upstream by sc_pileup (engine `popscle`) and
# shared with demuxlet; this stage no longer runs dsc-pileup or the per-chrom
# split/merge. freemuxlet resumes from the merged pileup prefix.

# freemuxlet clusters droplets into anonymous genotype groups de novo (no
# reference genotypes). Pure popscle, so it runs in the focused popscle image.
# Naming those clusters to donors is a separate Demuxafy R step (freemuxlet_assign)
# that needs R + the Demuxafy scripts, so it lives in its own rule/image.
rule freemuxlet:
    input:
        var=pileup_path("popscle", "pileup.var.gz"),
        cel=pileup_path("popscle", "pileup.cel.gz"),
        plp=pileup_path("popscle", "pileup.plp.gz"),
        umi=pileup_path("popscle", "pileup.umi.gz"),
        barcodes=lookup(captures, "barcodes")
    output:
        samples=out_dir/"freemuxlet/{capture_id}/freemuxlet.clust1.samples.gz",
        vcf=out_dir/"freemuxlet/{capture_id}/freemuxlet.clust1.vcf.gz",
        lmix=out_dir/"freemuxlet/{capture_id}/freemuxlet.lmix"
    params:
        N=get_number_of_samples_in_capture_pool
    container:
        config["containers"]["popscle"]
    log:
        logs/"freemuxlet/{capture_id}/output.freemuxlet.log"
    shell:
        """
        echo "Looking for {params.N} samples in capture"
        plp={input.plp}
        in_prefix=${{plp%.plp.gz}}
        lmix={output.lmix}
        out_prefix=${{lmix%.lmix}}
        popscle freemuxlet \
            --plp $in_prefix \
            --out $out_prefix \
            --group-list {input.barcodes} \
            --nsample {params.N} \
            > {log} 2>&1
        """


# Match freemuxlet's anonymous clusters to named donors by correlating the
# cluster genotypes against the known donor VCF (Demuxafy's Assign_Indiv_by_Geno.R).
rule freemuxlet_assign:
    input:
        vcf=rules.freemuxlet.output.vcf,
        known_genotype_vcf=rules.known_genotypes.output.vcf
    output:
        tsv=out_dir/"freemuxlet/{capture_id}/ref_clust_pearson_correlations.tsv",
        png=out_dir/"freemuxlet/{capture_id}/ref_clust_pearson_correlation.png",
        key=temp(out_dir/"freemuxlet/{capture_id}/Genotype_ID_key.txt")
    container:
        config["containers"]["demuxafy"]
    log:
        logs/"freemuxlet/{capture_id}/output.freemuxlet_assign.log"
    shell:
        """
        freemux_dir=$(dirname {input.vcf})
        Assign_Indiv_by_Geno.R \
            -r {input.known_genotype_vcf} \
            -c {input.vcf} \
            -o $freemux_dir \
            > {log} 2>&1
        """


rule freemuxlet_summary:
    input:
        samples=rules.freemuxlet.output.samples
    output:
        summary=out_dir/"freemuxlet/{capture_id}/summary.txt"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        bash Freemuxlet_summary.sh {input.samples} > {output.summary}
        """
