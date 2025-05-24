import os
import gzip
import pandas as pd
import logging

logging.basicConfig(filename=snakemake.log[0], filemode="w",
                    format="%(asctime)s [%(levelname)s] %(message)s", level=logging.INFO)

logging.info("Merging plp files")

sample_map_file = snakemake.input.sample_map
snp_map_file = snakemake.input.snp_map
output_file = snakemake.output.plp

# Temporary file to prevent Snakemake from using an incomplete output
temp_output_file = output_file + "_tmp.gz"

# Load sample map and SNP map
logging.info("Loading sample map")
full_sample_map = pd.read_csv(sample_map_file, sep='\t', compression='gzip', index_col=0)

logging.info("Loading SNP map")
snp_map = pd.read_csv(snp_map_file, sep='\t', compression='gzip')
print(snp_map.head())

first_file = True

for plp_file in snakemake.input.plp:
    chromosome = os.path.basename(plp_file).split(".")[0]
    logging.info("Processing plp file for %s", chromosome)

    with gzip.open(plp_file, "rt") as plp_f:
        df = (
            pd.read_csv(
                plp_f, 
                sep='\t',
                dtype = {
                    '#DROPLET_ID': 'Int64',
                    'SNP_ID': 'Int64',
                    'ALLELES': str,
                    'BASEQS': str
                }
            )
            .rename(columns={'SNP_ID': 'original_SNP_ID', '#DROPLET_ID': 'original_DROPLET_ID'})
        )
        

    # SNP Mapping
    snps = snp_map[snp_map['CHROM'] == chromosome].drop(columns=['CHROM']).reset_index()

    # Droplet Mapping
    if chromosome in full_sample_map.columns:
        droplets = full_sample_map[['#DROPLET_ID', chromosome]].rename(columns={chromosome: 'original_DROPLET_ID'})
        df = pd.merge(df, droplets, how="left", on="original_DROPLET_ID")

    df = pd.merge(df, snps, how="left", on="original_SNP_ID").rename(columns={'#SNP_ID': 'SNP_ID'})

    if first_file:
        df.to_csv(temp_output_file, sep='\t', index=False, compression='gzip')
        first_file = False
    else:
        prev_df = pd.read_csv(temp_output_file, sep='\t', compression='gzip')
        merged = pd.concat([prev_df, df], ignore_index=True)
        merged.to_csv(temp_output_file, sep='\t', index=False, compression='gzip')

# Safely move temp file to final location
os.rename(temp_output_file, output_file)

logging.info("Successfully saved final merged plp file")
