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

rule known_genotypes:
    input:
        vcf=config['deps']['vcf']
    output:
        vcf=out_dir/"souporcell/annotation/{capture_id}/known_genotypes.vcf"
    params:
        donor_ids=get_donors_in_capture
    log:
        logs/"souporcell/annotation/known_genotypes_{capture_id}.log"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        bcftools view {input.vcf} \
            --samples {params.donor_ids} \
            --force-samples \
            -Ov -o {output.vcf} \
            2> {log}
        """

rule extract_genotyped_donors:
    input:
        vcf=rules.known_genotypes.output.vcf
    output:
        genotyped_donors=out_dir/"demuxlet/{capture_id}/genotyped_donors.txt"
    log:
        logs/"souporcell/annotation/extract_genotyped_donors_{capture_id}.log"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        bcftools query -l {input.vcf} > {output.genotyped_donors}
        """

# rule get_all_genotyped_donors:
#     input:
#         vcf=config['deps']['vcf']
#     output:
#         full_donor_list=temp(out_dir/"all_genotyped_donors.txt")
#     container:
#         config['containers']['demuxafy']
#     shell:
#         """
#         bcftools query -l {input.vcf} > {output.full_donor_list}
#         """

# rule extract_donor_ids:
#     input:
#         samples=config['deps']['samples'] # used by get_donors_in_capture
#     output:
#         donors_in_capture=temp(out_dir/"demuxlet/{capture_id}/donors_in_capture.txt")
#     params:
#         donor_ids=get_donors_in_capture
#     shell:
#         """
#         echo {params.donor_ids} | tr ',' '\n' > {output.donors_in_capture}
#         """

# rule extract_genotyped_donors:
#     input:
#         full_donor_list=rules.get_all_genotyped_donors.output.full_donor_list,
#         donors_in_capture=rules.extract_donor_ids.output.donors_in_capture
#     output:
#         genotyped_donors=temp(out_dir/"demuxlet/{capture_id}/genotyped_donors.txt")
#     shell:
#         """
#             grep -Fx -f {input.donors_in_capture} {input.full_donor_list} > {output.genotyped_donors}
#         """
