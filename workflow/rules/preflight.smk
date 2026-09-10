# ==================================
# Pre-flight: per-capture donor present/absent screen
# ----------------------------------
# A clustering-INDEPENDENT coarse check of whether each declared pool donor is
# actually PRESENT in a capture, giving an at-a-glance "who is here / who dropped
# out" verdict the report surfaces at the very top of every page.
#
# It reads ONLY the sc_pileup cellSNP pileup and the pool's known genotypes -- NO
# demux dependency -- so it is a genuine pre-flight that runs early/in parallel. It
# is a QC FLAG only: nothing downstream prunes the panel or sets souporcell k from
# it (both tested and rejected).
#
# Deliberately UNGROUPED (not in the `postproc` job group): it fans out per capture
# rather than queueing behind per-capture post-processing.
# ==================================


rule present_absent_screen:
    """
    Per-donor coarse present/absent verdict for one capture. Private sites =
    biallelic panel sites where a donor is the sole alt-carrier among pool donors;
    per cell we count covered private sites with an alt read, and call the donor
    present when enough cells show several private alleles at once. See
    workflow/scripts/present_absent_screen.py.
    """
    input:
        known_genotypes_vcf=rules.known_genotypes.output.vcf,
        base_vcf=pileup_path("cellsnp", "cellSNP.base.vcf.gz"),
        ad_mtx=pileup_path("cellsnp", "cellSNP.tag.AD.mtx"),
        dp_mtx=pileup_path("cellsnp", "cellSNP.tag.DP.mtx"),
        barcodes=pileup_path("cellsnp", "cellSNP.samples.tsv"),
    output:
        # Per-capture preflight data lives in the DATA tree (out_dir), not the
        # reports dir -- reports/ is for published HTML + key aggregated summaries.
        tsv=out_dir/"preflight/{capture_id}.present_absent.tsv",
        # Per-barcode x genotyped-donor private-alt score (hit/cov counts). The
        # post-hoc negative-test diagnostic reads this to sanity-check the
        # barcodes a by-elimination-recovered (ungenotyped) donor was assigned.
        score_matrix=out_dir/"preflight/{capture_id}.private_score.tsv.gz",
    params:
        present_absent=config['params'].get('present_absent', {}),
        # Full declared roster (all pool donors incl. any ungenotyped one). The
        # screen compares this against the genotyped donors in known_genotypes.vcf
        # to (a) HARD-FAIL if >1 donor is ungenotyped (by-elimination handles only
        # one) and (b) log/flag the single ungenotyped donor as a
        # recover-by-elimination candidate that cannot be positively screened.
        donor_ids=get_donors_in_capture,
    conda:
        "../envs/embedding.yaml"
    log:
        logs/"preflight/{capture_id}.present_absent.log"
    script:
        "../scripts/present_absent_screen.py"


rule aggregate_present_absent:
    """
    Sentinel for the pre-flight present/absent branch. Mirrors the other
    aggregate_* rules over captures.index so `rule all` can pull it in.
    """
    input:
        expand(
            out_dir/"preflight/{capture_id}.present_absent.tsv",
            capture_id=captures.index.tolist()
        )
    output:
        checkpoint_file=out_dir/"preflight/all.done"
    shell:
        """
        echo {input} > {output.checkpoint_file}
        """


rule negative_test:
    """
    POST-HOC sanity check of a by-elimination donor recovery. For a pool with one
    UNGENOTYPED donor X (flagged by the present/absent screen), test whether the
    barcodes hier-assigned to X are a genuinely distinct donor or mislabeled
    known-donor / doublet cells -- using the per-barcode private-alt score matrix
    (own-vs-baseline hit-rate per known donor) plus per-barcode sex. Emits an n/a
    row for pools with no ungenotyped donor. See workflow/scripts/negative_test.py.
    """
    input:
        present_absent=rules.present_absent_screen.output.tsv,
        score_matrix=rules.present_absent_screen.output.score_matrix,
        hier=rules.hierarchical_consensus.output.hier,
        # RELATEDNESS axis: de-novo souporcell cluster consensus (cluster ids match
        # hier's `cluster` column) + the per-pool donor genotype panel. Used to tell
        # a mislabel (X==D) from a distinct relative of D at D's hom-alt sites.
        cluster_genotypes=rules.souporcell_consensus.output.cluster_genotypes,
        panel_vcf=panel_path("informative.vcf.gz"),
    output:
        # Per-capture data in the DATA tree; only the aggregated summary below is
        # published under reports/.
        tsv=out_dir/"preflight/{capture_id}.negative_test.tsv"
    conda:
        "../envs/embedding.yaml"
    log:
        logs/"preflight/{capture_id}.negative_test.log"
    script:
        "../scripts/negative_test.py"


rule aggregate_negative_test:
    """
    Concatenate the per-capture negative-test rows into one cohort table AND finalise
    the relatedness verdict across co-pooled captures. The per-capture step emits a
    provisional related/suspect call from its LOCAL het_share; here we pool by donor
    pair and rewrite each verdict from the best-powered (max) het_share, since
    souporcell's cluster consensus only collapses het->hom at low cell counts. See
    workflow/scripts/aggregate_negative_test.py.
    """
    input:
        expand(
            out_dir/"preflight/{capture_id}.negative_test.tsv",
            capture_id=captures.index.tolist()
        )
    output:
        # The ONE published artifact for this branch: a cohort summary in reports/.
        summary=report_dir/"preflight/negative_test_summary.tsv",
        checkpoint_file=out_dir/"preflight/negative_test.all.done"
    conda:
        "../envs/embedding.yaml"
    log:
        logs/"preflight/aggregate_negative_test.log"
    script:
        "../scripts/aggregate_negative_test.py"
