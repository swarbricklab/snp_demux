# popscle dsc-pileup is produced upstream by sc_pileup (engine `popscle`) and
# shared by demuxlet AND freemuxlet; this stage no longer runs dsc-pileup or the
# per-chromosome split/merge. demuxlet resumes from the merged pileup prefix.

rule demuxlet:
    input:
        var=pileup_path("popscle", "pileup.var.gz"),
        cel=pileup_path("popscle", "pileup.cel.gz"),
        plp=pileup_path("popscle", "pileup.plp.gz"),
        umi=pileup_path("popscle", "pileup.umi.gz"),
        vcf=rules.known_genotypes.output.vcf,
        barcodes=lookup(captures,"barcodes")
    output:
        best=out_dir/"demuxlet/{capture_id}/demuxlet.best"
    params:
        tag="GT",
        coeff=config['params']['demuxlet']['coeff'],
        offset=config['params']['demuxlet']['offset']
    container:
        config["containers"]["demuxafy"]
    log:
        logs/"demuxlet/{capture_id}/demuxlet.log"
    shell:
        """
        var={input.var}
        in_prefix=${{var%.var.gz}}
        best={output.best}
        out_prefix=${{best%.best}}
        popscle demuxlet \
            --plp $in_prefix \
            --vcf {input.vcf} \
            --field {params.tag} \
            --group-list {input.barcodes} \
            --geno-error-coeff {params.coeff} \
            --geno-error-offset {params.offset} \
            --out $out_prefix \
            > {log} 2>&1
        """

rule demuxlet_summary:
    input:
        demuxlet_best=rules.demuxlet.output.best
    output:
        summary=out_dir/"demuxlet/{capture_id}/summary.txt"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        bash Demuxlet_summary.sh {input.demuxlet_best} > {output.summary}
        """
