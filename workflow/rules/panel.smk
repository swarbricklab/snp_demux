# ==================================
# Demux-side imputation-quality filter on the donor panel
# ----------------------------------
# The upstream imputation is permissive (INFO/R2 carried through); the R2 cut
# that suits demultiplexing is owned here (params.imputation), where its effect
# on the calls can be measured. See workflow/scripts/filter_panel.sh.
# ==================================


rule filter_panel:
    input:
        vcf=panel_path("informative.vcf.gz"),
    output:
        vcf=temp(out_dir / "panel/{capture_id}/informative.filtered.vcf.gz"),
        idx=temp(out_dir / "panel/{capture_id}/informative.filtered.vcf.gz.csi"),
    log:
        logs / "panel/{capture_id}.filter.log",
    container:
        config["containers"]["bcftools"]
    params:
        min_r2=imputation_param("min_r2"),
        min_maf=imputation_param("min_maf"),
    script:
        "../scripts/filter_panel.sh"


# Opt-in benchmark (not part of `rule all`):
#   snakemake r2_sweep --workflow-profile profiles/workflow ...
# Per capture and threshold, how many donor-discriminating sites each donor keeps.
rule r2_sweep_capture:
    input:
        vcf=panel_path("informative.vcf.gz"),
    output:
        tsv=out_dir / "benchmark/r2_sweep/{capture_id}.tsv",
    log:
        logs / "benchmark/r2_sweep_{capture_id}.log",
    container:
        config["containers"]["vireosnp"]
    params:
        thresholds=imputation_param("sweep_r2"),
        donors=get_donors_in_capture,
    script:
        "../scripts/r2_sweep.py"


rule r2_sweep:
    input:
        tables=expand(
            out_dir / "benchmark/r2_sweep/{capture_id}.tsv",
            capture_id=captures.index.tolist(),
        ),
    output:
        tsv=out_dir / "benchmark/r2_sweep.tsv",
    log:
        logs / "benchmark/r2_sweep.log",
    container:
        config["containers"]["vireosnp"]
    script:
        "../scripts/r2_sweep.py"
