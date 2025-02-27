chromosomes=['chr1', 'chr2', 'chr3','chr4', 'chr5', 'chr6','chr7', 'chr8', 'chr9','chr10', 'chr11', 'chr12',
'chr13', 'chr14', 'chr15','chr16', 'chr17', 'chr18','chr19', 'chr20', 'chr21', 'chr22']

rule split_bam_by_chromosome:
    input:
        bam=lookup(captures, "bam")
    output:
        bam=temp(out_dir/"split/{capture_id}/bams/{chrom}.bam"),
        index=temp(out_dir/"split/{capture_id}/bams/{chrom}.bam.bai")
    log:
        logs/"split_bam/{capture_id}/{chrom}.log"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        samtools view -b {input.bam} {wildcards.chrom} > {output.bam} 2> {log}
        samtools index {output.bam} 2>> {log}
        """

rule split_vcf_by_chromosome:
    input:
        vcf=config['deps']['vcf']
    output:
        vcf=temp(out_dir/"split/{capture_id}/vcfs/{chrom}.vcf.gz"),
        index=temp(out_dir/"split/{capture_id}/vcfs/{chrom}.vcf.gz.csi")
    log:
        logs/"split_vcf/{capture_id}/{chrom}.log"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        bcftools view -r {wildcards.chrom} {input.vcf} -Oz -o {output.vcf} 2> {log}
        bcftools index {output.vcf} 2>> {log}
        """

rule get_all_genotyped_donors:
    input:
        vcf=config['deps']['vcf']
    output:
        full_donor_list=temp(out_dir/"all_genotyped_donors.txt")
    container:
        config['containers']['demuxafy']
    shell:
        """
        bcftools query -l {input.vcf} > {output.full_donor_list}
        """

rule extract_donor_ids:
    input:
        samples=config['deps']['samples'] # used by get_donors_in_capture
    output:
        donors_in_capture=temp(out_dir/"demuxlet/{capture_id}/donors_in_capture.txt")
    params:
        donor_ids=get_donors_in_capture
    shell:
        """
        echo {params.donor_ids} | tr ',' '\n' > {output.donors_in_capture}
        """

rule extract_genotyped_donors:
    input:
        full_donor_list=rules.get_all_genotyped_donors.output.full_donor_list,
        donors_in_capture=rules.extract_donor_ids.output.donors_in_capture
    output:
        genotyped_donors=temp(out_dir/"demuxlet/{capture_id}/genotyped_donors.txt")
    shell:
        """
            grep -Fx -f {input.donors_in_capture} {input.full_donor_list} > {output.genotyped_donors}
        """

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

rule merge_pileup:
    input:
        cel=expand(out_dir/"demuxlet/pileup/{capture_id}/{chrom}.cel.gz", chrom=chromosomes, capture_id="{capture_id}"),
        var=expand(out_dir/"demuxlet/pileup/{capture_id}/{chrom}.var.gz", chrom=chromosomes, capture_id="{capture_id}"),
        plp=expand(out_dir/"demuxlet/pileup/{capture_id}/{chrom}.plp.gz", chrom=chromosomes, capture_id="{capture_id}"),
        umi=expand(out_dir/"demuxlet/pileup/{capture_id}/{chrom}.umi.gz", chrom=chromosomes, capture_id="{capture_id}")
    output:
        cel=temp(out_dir/"demuxlet/{capture_id}/pileup.cel.gz"),
        var=temp(out_dir/"demuxlet/{capture_id}/pileup.var.gz"),
        plp=temp(out_dir/"demuxlet/{capture_id}/pileup.plp.gz"),
        umi=temp(out_dir/"demuxlet/{capture_id}/pileup.umi.gz")
    log:
        logs/"merge_pileup/{capture_id}/merge_cel.log"
    script:
        "../scripts/merge_pileup.py"

rule demuxlet:
    input:
        var=rules.merge_pileup.output.var,
        cel=rules.merge_pileup.output.cel,
        plp=rules.merge_pileup.output.plp,
        umi=rules.merge_pileup.output.umi,
        vcf=rules.subset_vcf.output.donor_subset,
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