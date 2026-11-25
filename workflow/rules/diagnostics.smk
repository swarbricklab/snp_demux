# ==================================
# Per-capture diagnostics (design doc §D reports, §E leave-one-out)
# ----------------------------------
# Additive QC producers that feed the capture report and act as trusted/flagged
# gates: donor distinctiveness + MAF-survival, variant-space embedding + per-
# cluster silhouette, and a score-level leave-one-out validator. None of these
# modify demuxafy.smk or the five existing aggregate targets; a new
# aggregate_diagnostics sentinel collects them for `rule all`.
# ==================================


rule report_distinctiveness:
    """
    Per-donor unique-site count (sole alt-carrier among pool donors) from the real
    genotypes, plus the MAF-survival curve of those unique sites in the cellSNP
    pileup. TSV + barplot PNG.
    """
    input:
        known_genotypes_vcf=rules.known_genotypes.output.vcf,
        ad_mtx=pileup_path("cellsnp", "cellSNP.tag.AD.mtx"),
        dp_mtx=pileup_path("cellsnp", "cellSNP.tag.DP.mtx"),
        base_vcf=pileup_path("cellsnp", "cellSNP.base.vcf.gz"),
    output:
        tsv=report_dir/"diagnostics/{capture_id}.distinctiveness.tsv",
        png=report_dir/"diagnostics/{capture_id}.distinctiveness.png",
    conda:
        "../envs/vireosnp.yaml"
    log:
        logs/"diagnostics/{capture_id}.distinctiveness.log"
    script:
        "../scripts/report_distinctiveness.py"


rule report_variant_embedding:
    """
    Cells x sites alt-fraction -> TruncatedSVD(30) -> UMAP, coloured by the de-novo
    souporcell cluster, with per-cluster + overall silhouette on the SVD space.
    """
    input:
        ad_mtx=pileup_path("cellsnp", "cellSNP.tag.AD.mtx"),
        dp_mtx=pileup_path("cellsnp", "cellSNP.tag.DP.mtx"),
        samples=pileup_path("cellsnp", "cellSNP.samples.tsv"),
        clusters=rules.souporcell_troublet.output.clusters,
    output:
        png=report_dir/"diagnostics/{capture_id}.embedding.png",
        tsv=report_dir/"diagnostics/{capture_id}.embedding.metrics.tsv",
    resources:
        # PBS memory used at most 4.0 GB over 17 captures; killed ~18x at 4 GB
        # while grouped (when it ran at the default). Start at 8 GB and double
        # per retry (retries: 2): 8, 16, 32 GB, so attempt 3 equals the old flat
        # 32 GB. A mem key for this rule in profiles/workflow set-resources would
        # override this on every attempt: keep it out.
        mem_mb=lambda wildcards, attempt: 8000 * 2 ** (attempt - 1),
    conda:
        "../envs/embedding.yaml"
    log:
        logs/"diagnostics/{capture_id}.embedding.log"
    script:
        "../scripts/report_variant_embedding.py"


rule leave_one_out:
    """
    Score-level leave-one-out over the vireosnp GPb_diff matrix: does each donor's
    best cluster survive as distinct when that donor is removed, or get absorbed?
    """
    input:
        scores_full=rules.vireosnp.output.scores_full,
    output:
        tsv=report_dir/"diagnostics/{capture_id}.leave_one_out.tsv",
    params:
        margin_min=config['params']['hierarchical']['margin_min'],
    conda:
        "../envs/vireosnp.yaml"
    log:
        logs/"diagnostics/{capture_id}.leave_one_out.log"
    script:
        "../scripts/leave_one_out.py"


rule aggregate_diagnostics:
    """
    Sentinel for the per-capture diagnostics branch (distinctiveness, embedding,
    leave-one-out). Mirrors the existing aggregate_* rules over captures.index.
    """
    input:
        expand(
            [
                report_dir/"diagnostics/{capture_id}.distinctiveness.tsv",
                report_dir/"diagnostics/{capture_id}.embedding.metrics.tsv",
                report_dir/"diagnostics/{capture_id}.leave_one_out.tsv",
            ],
            capture_id=captures.index.tolist()
        )
    output:
        checkpoint_file=out_dir/"diagnostics/all.done"
    shell:
        """
        echo {input} > {output.checkpoint_file}
        """
