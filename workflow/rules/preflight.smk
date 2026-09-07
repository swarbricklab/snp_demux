# ==================================
# Per-capture donor present/absent screen
# ----------------------------------
# A largely clustering-INDEPENDENT check of whether each declared pool donor is
# actually PRESENT in a capture, giving an at-a-glance "who is here / who dropped
# out" verdict that the report surfaces at the very top of every page.
#
# It reads the sc_pileup cellSNP pileup and the pool's known genotypes AND the
# per-cell hierarchical consensus: the demux assignment is used only to (a) call a
# confidently-assigned donor present outright and (b) estimate each donor's ambient
# background from cells NOT assigned to it, so a dominant donor's own reads cannot
# inflate its own null. That makes the screen demux-INFORMED rather than strictly
# pre-flight (intended) -- the private-site test still stands on its own for the
# rare donors the demux misses.
#
# Deliberately UNGROUPED (not in the `postproc` job group): it fans out per capture
# rather than queueing behind per-capture post-processing.
# ==================================


rule present_absent_screen:
    """
    Per-donor private-site present/absent verdict for one capture. Private sites =
    biallelic panel sites where a donor is the sole alt-carrier among pool donors;
    per cell we count covered private sites with an alt read and compare the tail
    against the ambient-only binomial null. See workflow/scripts/present_absent_screen.py.
    """
    input:
        known_genotypes_vcf=rules.known_genotypes.output.vcf,
        base_vcf=pileup_path("cellsnp", "cellSNP.base.vcf.gz"),
        ad_mtx=pileup_path("cellsnp", "cellSNP.tag.AD.mtx"),
        dp_mtx=pileup_path("cellsnp", "cellSNP.tag.DP.mtx"),
        # cellSNP barcodes give the AD/DP matrix COLUMN order, so the per-cell
        # demux assignment can be aligned onto the matrix columns.
        samples=pileup_path("cellsnp", "cellSNP.samples.tsv"),
        # Per-cell hierarchical consensus (barcode, donor, status, ...). Used to
        # estimate each donor's ambient background from cells NOT assigned to it,
        # so a dominant donor's own reads cannot inflate its own null. This makes
        # the screen demux-INFORMED rather than strictly pre-flight (intended).
        hier=rules.hierarchical_consensus.output.hier,
    output:
        tsv=report_dir/"preflight/{capture_id}.present_absent.tsv"
    params:
        present_absent=config['params'].get('present_absent', {}),
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
            report_dir/"preflight/{capture_id}.present_absent.tsv",
            capture_id=captures.index.tolist()
        )
    output:
        checkpoint_file=out_dir/"preflight/all.done"
    shell:
        """
        echo {input} > {output.checkpoint_file}
        """
