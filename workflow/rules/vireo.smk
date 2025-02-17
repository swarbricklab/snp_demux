
rule cellsnp_lite:
    input:
        vcf=config['deps']['vcf'],
        bam=lookup(captures,"bam"),
        barcodes=lookup(captures,"barcodes")
    output:
        base=temp(out_dir/"cellSNP/{capture_id}/cellSNP.base.vcf.gz"),
        samples=temp(out_dir/"cellSNP/{capture_id}/cellSNP.samples.tsv"),
        ad_mtx=temp(out_dir/"cellSNP/{capture_id}/cellSNP.tag.AD.mtx"),
        dp_mtx=temp(out_dir/"cellSNP/{capture_id}/cellSNP.tag.DP.mtx"),
        oth_mtx=temp(out_dir/"cellSNP/{capture_id}/cellSNP.tag.OTH.mtx")
    params:
        minMAF=config['params']['cellsnp_lite']['minMAF'],
        minCOUNT=config['params']['cellsnp_lite']['minCOUNT']
    log:
        logs/"cellsnp_lite/{capture_id}.log"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        outdir=$(dirname {output.base})
        cellsnp-lite --samFile {input.bam} \
            --regionsVCF {input.vcf} \
            --barcodeFile {input.barcodes} \
            --outDir $outdir \
            --nproc {resources.threads} \
            --minMAF 0.1 \
            --minCOUNT 20 \
            --gzip \
            > {log} 2>&1
        """

rule subset_vcf:
    input:
        vcf=config['deps']['vcf'],
        idx=config['deps']['vcf']+'.csi',
        base=rules.cellsnp_lite.output.base
    output:
        donor_subset=temp(out_dir/"vireo/{capture_id}/donor_subset.vcf")
    params:
        sample_ids=get_samples_in_capture
    log:
        logs/"subset_vcf/subset_{capture_id}.log"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        bcftools view {input.vcf} \
            --regions-file {input.base} \
            --samples {params.sample_ids} \
            --force-samples \
            -Ov -o {output.donor_subset} \
            2> {log}
        """

rule vireo:
    input:
        base=rules.cellsnp_lite.output.base,
        ad_mtx=rules.cellsnp_lite.output.ad_mtx,
        dp_mtx=rules.cellsnp_lite.output.dp_mtx,
        samples=rules.cellsnp_lite.output.samples,
        donor_vcf=rules.subset_vcf.output.donor_subset
    output:
        donor_ids=out_dir/"vireo/{capture_id}/donor_ids.tsv",
        fig=out_dir/"vireo/{capture_id}/fig_GT_distance_estimated.pdf",
        prob_doublet=out_dir/"vireo/{capture_id}/prob_doublet.tsv.gz",
        prob_singlet=out_dir/"vireo/{capture_id}/prob_singlet.tsv.gz",
        summary=out_dir/"vireo/{capture_id}/summary.tsv"
    params:
        FIELD="GT",
        vireo_out_dir=out_dir/"vireo/{capture_id}/"
    log:
        console=logs/"vireo/{capture_id}.console.log",
        vireo=logs/"vireo/{capture_id}_log.txt"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        base={input.base}
        cell_data=$(dirname $base)
        donor_ids={output.donor_ids}
        vireo_dir=$(dirname $donor_ids)
        cmd="vireo --cellData=$cell_data \
            --donorFile={input.donor_vcf} \
            --outDir=$vireo_dir \
            --genoTag={params.FIELD} \
            > {log.console} 2>&1"
        echo $cmd; eval $cmd
        mv $vireo_dir/_log.txt {log.vireo}
        """
