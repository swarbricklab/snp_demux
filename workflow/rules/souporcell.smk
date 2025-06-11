rule souporcell:
    input:
        bam=lookup(captures,"bam"),
        barcodes=lookup(captures,"barcodes"),
        genome_fasta=config["ref"]["souporcell"]["genome_fasta"],
        known_genotype_vcf=rules.known_genotypes.output.vcf
    output:
        fastqs=temp(out_dir/"souporcell/{capture_id}/fastqs.done"),
        clusters=out_dir/"souporcell/{capture_id}/clusters.tsv",
        cluster_genotypes=out_dir/"souporcell/{capture_id}/cluster_genotypes.vcf",
        ambient_rna=out_dir/"souporcell/{capture_id}/ambient_rna.txt",
        bam=temp(out_dir/"souporcell/{capture_id}/souporcell_minimap_tagged_sorted.bam"),
        bai=temp(out_dir/"souporcell/{capture_id}/souporcell_minimap_tagged_sorted.bam.bai"),
        alt=out_dir/"souporcell/{capture_id}/alt.mtx",
        ref=out_dir/"souporcell/{capture_id}/ref.mtx",
        bc=temp(out_dir/"souporcell/{capture_id}/barcodes.tsv"),
        bcf_err=temp(out_dir/"souporcell/{capture_id}/bcftools.err"),
        cl_done=temp(out_dir/"souporcell/{capture_id}/clustering.done"),
        cl_tmp=temp(out_dir/"souporcell/{capture_id}/clusters_tmp.tsv"),
        cl_err=temp(out_dir/"souporcell/{capture_id}/clusters.err"),
        consensus=temp(out_dir/"souporcell/{capture_id}/consensus.done"),
        doub_err=temp(out_dir/"souporcell/{capture_id}/doublets.err"),
        map_err=temp(out_dir/"souporcell/{capture_id}/minimap.err"),
        remap=temp(out_dir/"souporcell/{capture_id}/remapping.done"),
        retag=temp(out_dir/"souporcell/{capture_id}/retagging.done"),
        retag_err=temp(out_dir/"souporcell/{capture_id}/retag.err"),
        troublet=temp(out_dir/"souporcell/{capture_id}/troublet.done"),
        variants=temp(out_dir/"souporcell/{capture_id}/variants.done"),
        vartrix=temp(out_dir/"souporcell/{capture_id}/vartrix.done"),
        vcf=temp(out_dir/"souporcell/{capture_id}/souporcell_merged_sorted_vcf.vcf.gz"),
        tbi=temp(out_dir/"souporcell/{capture_id}/souporcell_merged_sorted_vcf.vcf.gz.tbi")
    params:
        sample_number=get_number_of_samples_in_capture_pool # number of samples expected in sc pool -> number of SNP clusters to generate. 
    container:
        "docker://cumulusprod/souporcell:2.5"
    log:
        logs/"souporcell/{capture_id}/output.souporcell.log"
    shell:
        """
        souporcell_dir="$(dirname {output.clusters})"
        number_of_known_genotypes=$(bcftools query -l {input.known_genotype_vcf} | wc -l)
        sample_names=$(bcftools query -l {input.known_genotype_vcf} | tr '\n' ' ')
        if [ "$number_of_known_genotypes" == {params.sample_number} ]; then
            echo "Running souporcell with known genotypes" > {log}
            cmd="souporcell_pipeline.py -i {input.bam} -b {input.barcodes} \
                -f {input.genome_fasta} -o $souporcell_dir \
                -t {resources.threads} -k {params.sample_number} \
                --known_genotypes {input.known_genotype_vcf} \
                --known_genotypes_sample_names $sample_names"
            echo $cmd >> {log}
            eval $cmd >> {log} 2>&1
            touch {output.bcf_err} {output.vcf} {output.tbi}
        else
            echo "Running souporcell without known genotypes"
            cmd="souporcell_pipeline.py -i {input.bam} -b {input.barcodes} \
                -f {input.genome_fasta} -o $souporcell_dir \
                -t {resources.threads} -k {params.sample_number}"
            echo $cmd >> {log}
            eval $cmd >> {log} 2>&1
        fi
        """
