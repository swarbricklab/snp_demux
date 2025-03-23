import os
import gzip
import pandas as pd
import logging

logging.basicConfig(filename=snakemake.log[0], filemode="w",
                    format="%(asctime)s [%(levelname)s] %(message)s", level=logging.INFO)

logging.info("Merging var files (in-memory)")

output_file = snakemake.output.var
snp_map_file = snakemake.output.snp_map

var_list = []

for var_file in snakemake.input.var:
    chromosome = os.path.basename(var_file).split(".")[0]
    logging.info(f"Processing var file for {chromosome}: {var_file}")

    with gzip.open(var_file, "rt") as var_f:
        df = pd.read_csv(
            var_f,
            sep='\t',
            dtype={
                '#SNP_ID': 'Int64',
                'CHROM': str,
                'POS': 'Int64',
                'REF': str,
                'ALT': str,
                'AF': 'Float64',
            }
        )
        df = df.rename(columns={'#SNP_ID': 'original_SNP_ID'})

    var_list.append(df)

# Merge all var entries
merged_var = pd.concat(var_list, ignore_index=True)
merged_var.index.name = '#SNP_ID'

# Extract SNP map after merge
snp_map = merged_var[['original_SNP_ID', 'CHROM']]
snp_map = snp_map.reset_index()

# Write var file with correct index
merged_var = merged_var.drop('original_SNP_ID', axis=1)
merged_var.to_csv(output_file, sep='\t', index=True, compression='gzip')

# Write SNP map
snp_map.to_csv(snp_map_file, sep='\t', index=False, compression='gzip')

logging.info("Successfully saved merged var file and SNP map")
