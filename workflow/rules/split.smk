# The per-chromosome BAM/VCF split rules were only needed to feed popscle
# dsc-pileup; that pileup now comes from sc_pileup, so they are gone. This file
# retains the donor-aware genotype helpers used by souporcell/demuxlet.


rule known_genotypes:
    input:
        # Per-pool informative panel from sc_pileup (donor GTs, already pool-subset
        # + informative-filtered). The --samples step below is now a consistency
        # guard / donor-ordering enforcement rather than the primary subset.
        vcf=panel_path("informative.vcf.gz"),
    output:
        vcf=out_dir / "souporcell/annotation/{capture_id}/known_genotypes.vcf",
    log:
        logs / "souporcell/annotation/known_genotypes_{capture_id}.log",
    container:
        config["containers"]["bcftools"]
    params:
        donor_ids=get_donors_in_capture,
    shell:
        """
        bcftools view {input.vcf} \
            --samples {params.donor_ids} \
            --force-samples \
            -Ov -o {output.vcf} \
            2>{log}
        """


rule extract_genotyped_donors:
    input:
        vcf=rules.known_genotypes.output.vcf,
    output:
        genotyped_donors=out_dir / "demuxlet/{capture_id}/genotyped_donors.txt",
    log:
        logs / "souporcell/annotation/extract_genotyped_donors_{capture_id}.log",
    container:
        config["containers"]["bcftools"]
    shell:
        """
        bcftools query -l {input.vcf} >{output.genotyped_donors}
        """
