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
            --minMAF {params.minMAF} \
            --minCOUNT {params.minCOUNT} \
            --gzip \
            > {log} 2>&1
        """

rule subset_vcf:
    input:
        vcf=config['deps']['vcf'],
        idx=config['deps']['vcf']+'.csi',
        base=rules.cellsnp_lite.output.base
    output:
        donor_subset=temp(out_dir/"vireo/tmp/{capture_id}/donor_subset.vcf")
    params:
        donor_ids=get_donors_in_capture
    log:
        logs/"subset_vcf/subset_{capture_id}.log"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        bcftools view {input.vcf} \
            --regions-file {input.base} \
            --samples {params.donor_ids} \
            --force-samples \
            -Ov -o {output.donor_subset} \
            2> {log}
        """

# This rule runs the vireo algorithm to estimate the probability of each cell being a doublet or singlet,
# and assigns donor IDs to cells based on their genotype data.

rule vireo:
    input:
        base=rules.cellsnp_lite.output.base,
        samples=rules.cellsnp_lite.output.samples,
        ad_mtx=rules.cellsnp_lite.output.ad_mtx,
        dp_mtx=rules.cellsnp_lite.output.dp_mtx,
        donor_vcf=rules.subset_vcf.output.donor_subset
    output:
        donor_ids=temp(out_dir/"vireo/tmp/{capture_id}/donor_ids.tsv"),
        fig_est=temp(out_dir/"vireo/tmp/{capture_id}/fig_GT_distance_estimated.pdf"),
        fig_input=temp(out_dir/"vireo/tmp/{capture_id}/fig_GT_distance_input.pdf"), # only for unknown genotypes
        prob_doublet=temp(out_dir/"vireo/tmp/{capture_id}/prob_doublet.tsv.gz"),
        prob_singlet=temp(out_dir/"vireo/tmp/{capture_id}/prob_singlet.tsv.gz"),
        prop_ambient=temp(out_dir/"vireo/tmp/{capture_id}/prop_ambient.tsv"),
        summary=temp(out_dir/"vireo/tmp/{capture_id}/summary.tsv"),
        vcf=temp(out_dir/"vireo/tmp/{capture_id}/GT_donors.vireo.vcf.gz") # only for unknown genotypes
    params:
        field="GT",
        n=get_number_of_samples_in_capture_pool
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
            --nDonor={params.n} \
            --nproc={resources.threads} \
            --callAmbientRNAs \
            --outDir=$vireo_dir \
            --genoTag={params.field} \
            > {log.console} 2>&1"
        echo $cmd; eval $cmd
        mv $vireo_dir/_log.txt {log.vireo}

        # Touch optional outputs if they don't exist
        [ -f {output.fig_input} ] || touch {output.fig_input}
        [ -f {output.vcf} ] || touch {output.vcf}
        [ -f {output.fig_est} ] || touch {output.fig_est}
        """

# This rule relabels the outputs of the vireo rule using donor information and writes them to the final output locations.
# It uses the relabel_vireo.R script, which takes all inputs and outputs from Snakemake variables.

rule relabel_vireo:
    input:
        samples_csv=config["deps"]["samples"],
        donor_ids=rules.vireo.output.donor_ids,
        prob_doublet=rules.vireo.output.prob_doublet,
        prob_singlet=rules.vireo.output.prob_singlet,
        prop_ambient=rules.vireo.output.prop_ambient,
        fig_est=rules.vireo.output.fig_est,
        fig_input=rules.vireo.output.fig_input,
        summary=rules.vireo.output.summary
    output:
        donor_ids=out_dir/"vireo/{capture_id}/donor_ids.tsv",
        sample_map=temp(out_dir/"vireo/{capture_id}/sample_map.tsv"),
        prob_doublet=out_dir/"vireo/{capture_id}/prob_doublet.tsv.gz",
        prob_singlet=out_dir/"vireo/{capture_id}/prob_singlet.tsv.gz",
        prop_ambient=out_dir/"vireo/{capture_id}/prop_ambient.tsv",
        fig_est_png=out_dir/"vireo/{capture_id}/fig_GT_distance_estimated.png",
        input_dist_png=out_dir/"vireo/{capture_id}/fig_GT_distance_input.png",
        summary=out_dir/"vireo/{capture_id}/summary.tsv"
    log:
        logs/"vireo/relabel_{capture_id}.log"
    conda:
        "../envs/relabel_vireo.yaml"
    script:
        "../scripts/relabel_vireo.R"

# This rule reheaders the VCF file produced by the vireo rule using the sample_map.tsv file from the relabel_vireo rule.
# It uses bcftools reheader to update the sample names in the VCF.

rule reheader_vireo_vcf:
    input:
        vcf=rules.vireo.output.vcf,
        sample_map=rules.relabel_vireo.output.sample_map
    output:
        vcf=out_dir/"vireo/{capture_id}/GT_donors.vireo.vcf.gz"
    log:
        logs/"vireo/reheader_{capture_id}.log"
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        if [ -s {input.vcf} ]; then
            bcftools reheader -s {input.sample_map} {input.vcf} -o {output.vcf} > {log} 2>&1
        else
            echo "No donor vcf file produced by vireo" > {output.vcf}
        fi
        """

# Summarise ambient RNA per barcode using vireo outputs
rule summarise_ambient:
    input:
        prop_ambient=rules.relabel_vireo.output.prop_ambient,
        donor_ids=rules.relabel_vireo.output.donor_ids
    output:
        summary=out_dir/"vireo/{capture_id}/ambient_summary.csv",
        stats=out_dir/"vireo/{capture_id}/ambient_stats.csv",
        plot=out_dir/"vireo/{capture_id}/ambient_plot.png"
    log:
        logs/"vireo/summarise_ambient_{capture_id}.log"
    conda:
        "../envs/ambient.yaml"
    script:
        "../scripts/summarise_ambient.R"
