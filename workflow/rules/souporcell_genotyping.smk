# Use vireoSNP to do the scRNA SNP to known genotype mapping
rule vireosnp:
    input:
        clusters=out_dir/"souporcell/{capture_id}/clusters.tsv",
        cluster_genotypes=out_dir/"souporcell/{capture_id}/cluster_genotypes.vcf",
        genotyped_donors=rules.extract_genotyped_donors.output.genotyped_donors,
        known_genotype_vcf=rules.known_genotypes.output.vcf
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
        config["containers"]["vireosnp"]
    log:
        logs/"vireosnp/{capture_id}.log"
    script:
        "../scripts/annotate.py"
