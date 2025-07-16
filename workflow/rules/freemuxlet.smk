rule dsc_pileup_freemuxlet_by_chromosome:
    input:
        vcf=rules.split_vcf_by_chromosome.output.vcf,
        bam=rules.split_bam_by_chromosome.output.bam,
        barcodes=lookup(captures, "barcodes")
    output:
        cel=temp(out_dir/"freemuxlet/pileup/{capture_id}/{chrom}.cel.gz"),
        plp=temp(out_dir/"freemuxlet/pileup/{capture_id}/{chrom}.plp.gz"),
        umi=temp(out_dir/"freemuxlet/pileup/{capture_id}/{chrom}.umi.gz"),
        var=temp(out_dir/"freemuxlet/pileup/{capture_id}/{chrom}.var.gz")
    container:
        config["containers"]["demuxafy"]
    log:
        logs/"dsc_pileup_freemuxlet/{capture_id}/{chrom}.log"
    shell:
        """
        cel={output.cel}
        out_prefix=${{cel%.cel.gz}}
        popscle dsc-pileup \
            --sam {input.bam} \
            --vcf {input.vcf} \
            --group-list {input.barcodes} \
            --out $out_prefix \
            > {log} 2>&1
        """

rule merge_cel_freemuxlet:
    input:
        cel=expand(out_dir/"freemuxlet/pileup/{capture_id}/{chrom}.cel.gz", chrom=chromosomes, capture_id="{capture_id}")
    output:
        cel=temp(out_dir/"freemuxlet/{capture_id}/pileup.cel.gz"),
        sample_map=temp(out_dir/"freemuxlet/{capture_id}/full_sample_map.tsv.gz")
    log:
        logs/"merge_pileup/freemuxlet/{capture_id}/merge_cel.log"
    script:
        "../scripts/merge_cel.py"

rule merge_var_freemuxlet:
    input:
        var=expand(out_dir/"freemuxlet/pileup/{capture_id}/{chrom}.var.gz", chrom=chromosomes, capture_id="{capture_id}")
    output:
        var=temp(out_dir/"freemuxlet/{capture_id}/pileup.var.gz"),
        snp_map=temp(out_dir/"freemuxlet/{capture_id}/snp_map.tsv.gz")
    log:
        logs/"merge_pileup/freemuxlet/{capture_id}/merge_var.log"
    script:
        "../scripts/merge_var.py"


rule wrangle_plp_freemuxlet:
    input:
        plp=out_dir/"freemuxlet/pileup/{capture_id}/{chrom}.plp.gz",
        sample_map=rules.merge_cel_freemuxlet.output.sample_map,
        snp_map=rules.merge_var_freemuxlet.output.snp_map
    output:
        plp=temp(out_dir/"freemuxlet/wrangle/{capture_id}/{chrom}.plp")
    log:
        logs/"merge_pileup/freemuxlet/{capture_id}/wrangle_plp_{chrom}.log"
    script:
        "../scripts/wrangle_plp.py"

rule merge_plp_freemuxlet:
    input:
        plp=expand(out_dir/"freemuxlet/wrangle/{capture_id}/{chrom}.plp", 
                chrom=chromosomes, 
                capture_id="{capture_id}"
            )
    output:
        plp=temp(out_dir/"freemuxlet/{capture_id}/pileup.plp.gz")
    log:
        logs/"merge_pileup/freemuxlet/{capture_id}/merge_plp.log"
    shell:
        """
        csvstack --tabs {input.plp} | csvformat -T | gzip -n > {output.plp}
        """

rule wrangle_umi_freemuxlet:
    input:
        umi=out_dir/"freemuxlet/pileup/{capture_id}/{chrom}.umi.gz",
        sample_map=rules.merge_cel_freemuxlet.output.sample_map
    output:
        umi=temp(out_dir/"freemuxlet/wrangle/{capture_id}/{chrom}.umi.gz")
    log:
        logs/"merge_pileup/freemuxlet/{capture_id}/wrangle_umi_{chrom}.log"
    script:
        "../scripts/wrangle_umi.py"

rule merge_umi_freemuxlet:
    input:
        umi=expand(out_dir/"freemuxlet/wrangle/{capture_id}/{chrom}.umi.gz", chrom=chromosomes, capture_id="{capture_id}")
    output:
        umi=temp(out_dir/"freemuxlet/{capture_id}/pileup.umi.gz")
    log:
        logs/"merge_pileup/freemuxlet/{capture_id}/merge_umi.log"
    shell:
        "zcat {input.umi} | gzip > {output.umi}"

rule freemuxlet:
    input:
        var=rules.merge_var_freemuxlet.output.var,
        cel=rules.merge_cel_freemuxlet.output.cel,
        plp=rules.merge_plp_freemuxlet.output.plp,
        umi=rules.merge_umi_freemuxlet.output.umi,
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
