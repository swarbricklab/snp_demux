# ==================================
# freemuxlet genotyping
# ----------------------------------
# Run to Correlate Cluster to Donor Reference SNP Genotypes
#
# TODO: This is not part of the DAG, where does it interact with other rules? 
# ==================================
rule freemuxlet_genotyping:
    input:
        VCF=imputation/"demuxafy_input/reheader_Merged_Imputed_R2_0.3_MAF0.05_Exon_sorted.vcf"
    output:
        Genotype_ID_key=freemuxlet/"{capture_id}/Genotype_ID_key.txt"
    params:
        freemuxlet_out_dir=freemuxlet/"{capture_id}/"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        Assign_Indiv_by_Geno.R -r {input.VCF} -c {params.freemuxlet_out_dir}freemuxlet.clust1.vcf.gz -o {params.freemuxlet_out_dir}
        """