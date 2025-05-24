import os
import gzip
import pandas as pd
import logging

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    format="%(asctime)s [%(levelname)s] %(message)s",
    level=logging.INFO,
)

chromosome = snakemake.wildcards.chrom
logging.info(f"Wrangling PLP for {chromosome}")

sample_map_file = snakemake.input.sample_map
snp_map_file = snakemake.input.snp_map
input_file = snakemake.input.plp
output_file = snakemake.output.plp

# Load maps
sample_map = pd.read_csv(sample_map_file, sep='\t', compression='gzip', index_col=0)
snp_map = pd.read_csv(snp_map_file, sep='\t', compression='gzip')

# Load and rename columns in plp
with gzip.open(input_file, "rt") as plp_file:
    df = (
            pd.read_csv(
                plp_file, 
                sep='\t',
                dtype = {
                    '#DROPLET_ID': 'Int64',
                    'SNP_ID': 'Int64',
                    'ALLELES': str,
                    'BASEQS': str
                }
            )
            .rename(
                columns={
                    'SNP_ID': 'original_SNP_ID', 
                    '#DROPLET_ID': 'original_DROPLET_ID'
                    }
            )
    )

print(df.head())

# Extract chromosome-specific SNPs
snps = (
    snp_map[snp_map['CHROM'] == chromosome]
    .drop(columns=['CHROM'])
)
print(snps.head())

droplets = (
    sample_map[['#DROPLET_ID', chromosome]]
    .rename(columns={chromosome: 'original_DROPLET_ID'})
)
print(droplets.head())

# Perform join and reduce to required columns
plp_joined = (
    pd.merge(df, snps, how='left', on='original_SNP_ID')
      .merge(droplets, how='left', on='original_DROPLET_ID')
      .filter(items=['#DROPLET_ID', '#SNP_ID', 'ALLELES', 'BASEQS'])
      .rename(columns={'#SNP_ID': 'SNP_ID'})
)
print(plp_joined.head())

# Save output
plp_joined.to_csv(output_file, sep='\t', index=False)
logging.info(f"Wrote wrangled PLP for {chromosome} to {output_file}")
