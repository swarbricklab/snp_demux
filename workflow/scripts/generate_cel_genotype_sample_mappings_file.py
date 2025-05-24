# ===================
# title: "generate_cel_genotype_sample_mappings_file"
# description: Generates Genotype-sample ID mapping list file from the CEL files in a VCF header
# project: "workflow-demuxafy"
# author: "Daniel Roden"
# email: "d.roden@garvan.org.au"
# company: "Garvan Institute"
# created: 07/06/2023
# -------------------
# NOTE: this is a snakmake script that uses the "snakemake" object internals
# ===================

import logging as log
import pysam
from pysam import VariantFile
import pandas as pd
import os

# extract just the CEL file-name from the full-path
def get_cel_filename(cel_filename):
    if pd.notna(cel_filename):
        cel_filename = os.path.basename(cel_filename)
    return cel_filename


log.basicConfig(filename = snakemake.log[0], 
                    level = log.INFO, 
                    format = '%(asctime)s:%(levelname)s:%(message)s')                    
log.info("LOG FILE: " + snakemake.log[0])

log.info("----------------------------------") 
log.info("CONFIG")

log.info("INPUT")
file_samples = snakemake.input["samples"]
log.info("file_samples: " + file_samples)

file_vcf = snakemake.input["vcf"]
log.info("file_vcf: " + file_vcf)
log.info("")

log.info("OUTPUT")
file_genotype_sample_mappings = snakemake.output["genotype_sample_mappings"]
log.info("file_genotype_sample_mappings: " + file_genotype_sample_mappings)
log.info("")

log.info("WILDCARDS")
log.info("")

log.info("PARAMS")
log.info("")

log.info("----------------------------------") 

log.info("get the samples ")
df_samples = pd.read_csv(file_samples, sep="\t", dtype=str)

log.info("extract just the CEL file-name from the full-path")
df_samples["cel_filename"] = df_samples["genotype_cel_file"].apply(get_cel_filename)

log.info("get the VCF samples")
vcf_in = VariantFile(file_vcf)
sample_list = vcf_in.header.samples

sample_list = [str(i) for i in sample_list]
#sample_list = sample_list.tolist()
log.info("sample_list")
log.info(sample_list)

log.info("Re-map existing VCF CEL file names to genotype_sample_id")

df_sample_mapping = df_samples.loc[df_samples["cel_filename"].isin(sample_list), ["cel_filename", "GENOTYPE_SAMPLE_ID"]]
df_sample_mapping = df_sample_mapping.drop_duplicates()

log.info("Mapped CEL - sample_list")
log.info(df_sample_mapping)

# A more robust version may be to use the old -> new mapping format
# see: https://samtools.github.io/bcftools/bcftools.html#reheader
log.info("write re-mapped IDs to file")
df_sample_mapping.to_csv(file_genotype_sample_mappings, sep = " ", index = False)

log.info("")
log.info("PROCESSING COMPLETE")
