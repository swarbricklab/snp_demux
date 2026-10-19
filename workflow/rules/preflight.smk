# Read-only preflight: validate each capture's informative donor panel against
# the input requirements in docs/data_prep_requirements.md, so a misprepared
# panel is caught up front instead of as a cryptic popscle/demuxlet FATAL.
#
# NOT part of `rule all` or the demux DAG — it reads nothing the demux rules
# write and writes nothing they read, so running (or skipping) it cannot change
# any demux result. Run it explicitly before a full run:
#
#   snakemake check_inputs --workflow-profile profiles/workflow [--use-singularity ...]


rule check_panel:
    input:
        vcf=panel_path("informative.vcf.gz"),
    output:
        report=out_dir / "preflight/{capture_id}.panel_check.txt",
    log:
        logs / "preflight/{capture_id}.log",
    container:
        config["containers"]["bcftools"]
    script:
        "../scripts/check_panel.sh"


rule check_inputs:
    input:
        expand(
            out_dir / "preflight/{capture_id}.panel_check.txt",
            capture_id=captures.index.tolist(),
        ),
    output:
        out_dir / "preflight/ALL_PANELS_OK.txt",
    log:
        logs / "preflight/aggregate.log",
    container:
        config["containers"]["bcftools"]
    shell:
        "cat {input} > {output} 2> {log} && echo '# All panels passed preflight.' >> {output}"
