# The per-chromosome BAM/VCF split rules were only needed to feed popscle
# dsc-pileup; that pileup now comes from sc_pileup, so they are gone. This file
# retains the donor-aware genotype helpers used by souporcell/demuxlet.

rule known_genotypes:
    input:
        # Per-pool informative panel from sc_pileup (donor GTs, already pool-subset
        # + informative-filtered). A donor absent from every genotype source is
        # OMITTED from this panel by build_harmonized_panel.sh (NOT an all-./.
        # column). The --samples step below selects only the pool donors that are
        # ACTUALLY PRESENT as samples in the panel (declared donors ∩ panel
        # samples), preserving declared order. It deliberately does NOT use
        # --force-samples: that would RE-ADD the omitted donor as an all-./. column,
        # making it look "genotyped" to extract_genotyped_donors and undoing the
        # by-elimination recovery + causing vireo cluster mis-steal. Mirrors the
        # name-filter guard in souporcell.smk. The omitted donor still appears in
        # the samples.csv roster annotate.py diffs against, so it becomes the "one
        # unknown" and is recovered downstream.
        vcf=panel_path("informative.vcf.gz")
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
        # Intersect declared pool donors with the panel's actual samples (order
        # preserved). strict-mode-safe: `|| true` swallows any grep SIGPIPE; the
        # if/then keeps a not-found grep from aborting the loop under `set -e`.
        avail=$(bcftools query -l {input.vcf}) || true
        keep=""
        for d in $(echo {params.donor_ids} | tr ',' ' '); do
            if printf '%s\\n' "$avail" | grep -qxF -- "$d"; then keep="$keep,$d"; fi
        done
        keep=${{keep#,}}
        echo "known_genotypes {wildcards.capture_id}: declared=[{params.donor_ids}] kept=[$keep]" > {log}
        bcftools view {input.vcf} \
            --samples "$keep" \
            -Ov -o {output.vcf} \
            2>> {log}
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
