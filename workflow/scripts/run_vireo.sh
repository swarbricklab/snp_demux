#! /bin/bash

# This script runs vireo 
# As well as wraps the sample list - which is an input for vireo 
set -e

# INPUT 
sample_list=snakemake.input["samples"]
#log=${snakemake.log[0]}
#bam=snakemake.input["bam"]
#barcodes=snakemake.input["barcodes"]
#vireo_out_dir=snakemake.params["vireo_out_dir"]
#VCF=snakemake.input["VCF"]
#vcf_zip=snakemake.input["vcf_zip"]
#FIELD=snakemake.params["FIELD"]

#sample_list=$(tr -s '\n ' ',' < sample_list) 
#sample_list=${sample_list::-1}
echo "sample list: $sample_list"
#echo $bam
#echo $barcodes 
#echo $vireo_out_dir

#touch $log 
#echo "sample list: $sample_list" >  $log

#sample_list > "$vireo_out_dir/test.txt"
# OUTPUT 


#cellsnp-lite -s {input.bam} -b {input.barcodes} -O {params.vireo_out_dir} -R {input.VCF} -p 20 --minMAF 0.1 --minCOUNT 20 --gzip
#bcftools view {input.vcf_zip} -R {params.vireo_out_dir}cellSNP.base.vcf.gz -s $sample_list -Ov -o {params.vireo_out_dir}donor_subset.vcf
#vireo -c {params.vireo_out_dir} -d {params.vireo_out_dir}donor_subset.vcf -o {params.vireo_out_dir} -t {params.FIELD}