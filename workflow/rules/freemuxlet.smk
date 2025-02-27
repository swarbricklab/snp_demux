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

rule merge_pileup_freemuxlet:
    input:
        cel=expand(out_dir/"freemuxlet/pileup/{capture_id}/{chrom}.cel.gz", chrom=chromosomes, capture_id="{capture_id}"),
        var=expand(out_dir/"freemuxlet/pileup/{capture_id}/{chrom}.var.gz", chrom=chromosomes, capture_id="{capture_id}"),
        plp=expand(out_dir/"freemuxlet/pileup/{capture_id}/{chrom}.plp.gz", chrom=chromosomes, capture_id="{capture_id}"),
        umi=expand(out_dir/"freemuxlet/pileup/{capture_id}/{chrom}.umi.gz", chrom=chromosomes, capture_id="{capture_id}")
    output:
        cel=temp(out_dir/"freemuxlet/{capture_id}/pileup.cel.gz"),
        var=temp(out_dir/"freemuxlet/{capture_id}/pileup.var.gz"),
        plp=temp(out_dir/"freemuxlet/{capture_id}/pileup.plp.gz"),
        umi=temp(out_dir/"freemuxlet/{capture_id}/pileup.umi.gz")
    log:
        logs/"merge_pileup_freemuxlet/{capture_id}/merge_cel.log"
    script:
        "../scripts/merge_pileup.py"


rule freemuxlet:
    input:
        var=rules.merge_pileup_freemuxlet.output.var,
        cel=rules.merge_pileup_freemuxlet.output.cel,
        plp=rules.merge_pileup_freemuxlet.output.plp,
        barcodes=lookup(captures, "barcodes")
    output:
        samples=out_dir/"freemuxlet/{capture_id}/freemuxlet.clust1.samples.gz",
        vcf=out_dir/"freemuxlet/{capture_id}/freemuxlet.clust1.vcf.gz",
        lmix=out_dir/"freemuxlet/{capture_id}/freemuxlet.lmix"
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
