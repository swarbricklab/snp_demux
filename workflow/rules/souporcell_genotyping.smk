# ==================================
# souporcell genotyping
# ----------------------------------
# Run to Correlate Cluster to Donor Reference SNP Genotypes
# ==================================
# TODO: is this redundant if we use the demuxafy::combined_results.R "--souporcell" and asssociated "--souporcell_correlation_limit" procssing
rule souporcell_genotyping:
    input:
        VCF=config['deps']['vcf'],
        souporcell=out_dir/"souporcell/{capture_id}/cluster_genotypes.vcf"
    output:
        Genotype_ID_key=out_dir/"souporcell/{capture_id}/Genotype_ID_key.txt"
    container:
        config["containers"]["demuxafy"]
    log:
        logs/"souporcell_genotyping/{capture_id}.log"
    shell:
        """
        souporcell_dir=$(dirname {output.Genotype_ID_key})
        Assign_Indiv_by_Geno.R -r {input.VCF} -c {input.souporcell} -o $souporcell_dir > {log} 2>&1
        """

# ==================================
# run vireosnp
# ----------------------------------
# Use vireoSNP to do the scRNA SNP to known genotype mapping
# NOTE: this was taken from here: https://git.gimr.garvan.org.au/CTP/soup-or-cell/-/blob/main/workflow/rules/annotation.smk
# ==================================
rule vireosnp:
    input:
        clusters=out_dir/"souporcell/{capture_id}/clusters.tsv",
        cluster_genotypes=out_dir/"souporcell/{capture_id}/cluster_genotypes.vcf",
        # NOTE: this has to be int encoded otherwise it generates an error
        known_genotype_vcf=config['deps']['vcf']
    output:
        clusters_annotated=out_dir/"souporcell/annotation/vireosnp/{capture_id}/clusters_annotated.tsv",
        scores_matched=out_dir/"souporcell/annotation/vireosnp/{capture_id}/matched_GPb_diff_scores_cluster_genotype.tsv",
        scores_full=out_dir/"souporcell/annotation/vireosnp/{capture_id}/full_GPb_diff_scores_cluster_genotype.tsv",
        matching_heatmap_pdf=report(out_dir/"souporcell/annotation/vireosnp/{capture_id}/matched_GPb_diff_scores_cluster_genotype.heatmap.pdf", caption="../report/vireosnp.rst", category="{projectname} ({capture_id})", subcategory="Matching pairs"),
        matching_heatmap_png=out_dir/"souporcell/annotation/vireosnp/{capture_id}/matched_GPb_diff_scores_cluster_genotype.heatmap.png",
        full_heatmap_pdf=report(out_dir/"souporcell/annotation/vireosnp/{capture_id}/full_GPb_diff_scores_cluster_genotype.heatmap.pdf", caption="../report/vireosnp.rst", category="({capture_id})", subcategory="All pairs"),
        full_heatmap_png=out_dir/"souporcell/annotation/vireosnp/{capture_id}/full_GPb_diff_scores_cluster_genotype.heatmap.png",
        grouped_scores_boxplot=out_dir/"souporcell/annotation/vireosnp/{capture_id}/grouped_GPb_diff_scores.boxplot.pdf",
        final_sample_cluster_mappings=out_dir/"souporcell/annotation/vireosnp/{capture_id}/final_sample_to_snpcluster_mappings.tsv"
    params:
        dir_output=out_dir/"souporcell/annotation/vireosnp/{capture_id}/"
    container:
        "docker://swarbricklab/vireo_snp:0.5.6"
    log:
        logs/"vireosnp/{capture_id}.log"
    script:
        "../scripts/annotate.py"
