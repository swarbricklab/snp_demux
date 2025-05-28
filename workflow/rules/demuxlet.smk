rule dsc_pileup_demuxlet_by_chromosome:
    input:
        bam=rules.split_bam_by_chromosome.output.bam,
        vcf=rules.split_vcf_by_chromosome.output.vcf,
        barcodes=lookup(captures, "barcodes"),
        donor_ids=rules.extract_genotyped_donors.output.genotyped_donors
    output:
        cel=temp(out_dir/"demuxlet/pileup/{capture_id}/{chrom}.cel.gz"),
        plp=temp(out_dir/"demuxlet/pileup/{capture_id}/{chrom}.plp.gz"),
        umi=temp(out_dir/"demuxlet/pileup/{capture_id}/{chrom}.umi.gz"),
        var=temp(out_dir/"demuxlet/pileup/{capture_id}/{chrom}.var.gz")
    log:
        logs/"dsc_pileup/{capture_id}/{chrom}.log"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        cel={output.cel}
        out_prefix=${{cel%.cel.gz}}
        popscle dsc-pileup \
            --sam {input.bam} \
            --vcf {input.vcf} \
            --group-list {input.barcodes} \
            --out $out_prefix \
            --sm-list {input.donor_ids} \
            > {log} 2>&1
        """

rule merge_cel_demuxlet:
    input:
        cel=expand(out_dir/"demuxlet/pileup/{capture_id}/{chrom}.cel.gz", chrom=chromosomes, capture_id="{capture_id}")
    output:
        cel=temp(out_dir/"demuxlet/{capture_id}/pileup.cel.gz"),
        sample_map=temp(out_dir/"demuxlet/{capture_id}/full_sample_map.tsv.gz")
    log:
        logs/"merge_pileup/demuxlet/{capture_id}/merge_cel.log"
    script:
        "../scripts/merge_cel.py"

rule merge_var_demuxlet:
    input:
        var=expand(out_dir/"demuxlet/pileup/{capture_id}/{chrom}.var.gz", chrom=chromosomes, capture_id="{capture_id}")
    output:
        var=temp(out_dir/"demuxlet/{capture_id}/pileup.var.gz"),
        snp_map=temp(out_dir/"demuxlet/{capture_id}/snp_map.tsv.gz")
    log:
        logs/"merge_pileup/demuxlet/{capture_id}/merge_var.log"
    script:
        "../scripts/merge_var.py"

rule wrangle_plp_demuxlet:
    input:
        plp=out_dir/"demuxlet/pileup/{capture_id}/{chrom}.plp.gz",
        sample_map=rules.merge_cel_demuxlet.output.sample_map,
        snp_map=rules.merge_var_demuxlet.output.snp_map
    output:
        plp=temp(out_dir/"demuxlet/wrangle/{capture_id}/{chrom}.plp")
    log:
        logs/"merge_pileup/demuxlet/{capture_id}/wrangle_plp_{chrom}.log"
    script:
        "../scripts/wrangle_plp.py"

rule merge_plp_demuxlet:
    input:
        plp=expand(out_dir/"demuxlet/wrangle/{capture_id}/{chrom}.plp", 
                chrom=chromosomes, 
                capture_id="{capture_id}"
            )
    output:
        plp=temp(out_dir/"demuxlet/{capture_id}/pileup.plp.gz")
    container:
        "docker://johnyaku/csvkit:1.1.1"
    log:
        logs/"merge_pileup/demuxlet/{capture_id}/merge_plp.log"
    shell:
        """
        csvstack --tabs {input.plp} | csvformat -T | gzip -n > {output.plp}
        """

rule wrangle_umi_demuxlet:
    input:
        umi=out_dir/"demuxlet/pileup/{capture_id}/{chrom}.umi.gz",
        sample_map=rules.merge_cel_demuxlet.output.sample_map
    output:
        umi=temp(out_dir/"demuxlet/wrangle/{capture_id}/{chrom}.umi.gz")
    log:
        logs/"merge_pileup/demuxlet/{capture_id}/wrangle_umi_{chrom}.log"
    script:
        "../scripts/wrangle_umi.py"

rule merge_umi_demuxlet:
    input:
        umi=expand(out_dir/"demuxlet/wrangle/{capture_id}/{chrom}.umi.gz", chrom=chromosomes, capture_id="{capture_id}")
    output:
        umi=temp(out_dir/"demuxlet/{capture_id}/pileup.umi.gz")
    log:
        logs/"merge_pileup/demuxlet/{capture_id}/merge_umi.log"
    shell:
        "zcat {input.umi} | gzip > {output.umi}"


rule demuxlet:
    input:
        var=rules.merge_var_demuxlet.output.var,
        cel=rules.merge_cel_demuxlet.output.cel,
        plp=rules.merge_plp_demuxlet.output.plp,
        umi=rules.merge_umi_demuxlet.output.umi,
        vcf=rules.known_genotypes.output.vcf,
        barcodes=lookup(captures, "barcodes")
    output:
        best=out_dir/"demuxlet/{capture_id}/demuxlet.best"
    params:
        tag="GT",
        coeff=config['params']['demuxlet']['coeff'],
        offset=config['params']['demuxlet']['offset']
    container:
        config["containers"]["demuxafy"]
    log:
        logs/"demuxlet/{capture_id}/demuxlet.log"
    shell:
        """
        var={input.var}
        in_prefix=${{var%.var.gz}}
        best={output.best}
        out_prefix=${{best%.best}}
        popscle demuxlet \
            --plp $in_prefix \
            --vcf {input.vcf} \
            --field {params.tag} \
            --group-list {input.barcodes} \
            --geno-error-coeff {params.coeff} \
            --geno-error-offset {params.offset} \
            --out $out_prefix \
            > {log} 2>&1
        """

rule demuxlet_summary:
    input:
        demuxlet_best=rules.demuxlet.output.best
    output:
        summary=out_dir/"demuxlet/{capture_id}/summary.txt"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        bash Demuxlet_summary.sh {input.demuxlet_best} > {output.summary}
        """