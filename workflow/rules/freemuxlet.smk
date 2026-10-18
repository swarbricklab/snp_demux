# popscle dsc-pileup is produced upstream by sc_pileup (engine `popscle`) and
# shared with demuxlet; this stage no longer runs dsc-pileup or the per-chrom
# split/merge. freemuxlet resumes from the merged pileup prefix.

rule freemuxlet:
    input:
        var=pileup_path("popscle", "pileup.var.gz"),
        cel=pileup_path("popscle", "pileup.cel.gz"),
        plp=pileup_path("popscle", "pileup.plp.gz"),
        umi=pileup_path("popscle", "pileup.umi.gz"),
        known_genotype_vcf=rules.known_genotypes.output.vcf,
        barcodes=lookup(captures, "barcodes")
    output:
        samples=out_dir/"freemuxlet/{capture_id}/freemuxlet.clust1.samples.gz",
        vcf=out_dir/"freemuxlet/{capture_id}/freemuxlet.clust1.vcf.gz",
        lmix=out_dir/"freemuxlet/{capture_id}/freemuxlet.lmix",
        tsv=out_dir/"freemuxlet/{capture_id}/ref_clust_pearson_correlations.tsv",
        png=out_dir/"freemuxlet/{capture_id}/ref_clust_pearson_correlation.png",
        key=temp(out_dir/"freemuxlet/{capture_id}/Genotype_ID_key.txt")
    params:
        N=get_number_of_samples_in_capture_pool
    container:
        config["containers"]["demuxafy"]
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

        freemux_dir=$(dirname {output.vcf})
        Assign_Indiv_by_Geno.R \
            -r {input.known_genotype_vcf} \
            -c {output.vcf} \
            -o $freemux_dir \
            >> {log} 2>&1
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
