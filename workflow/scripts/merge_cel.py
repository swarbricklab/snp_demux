import os
import gzip
import pandas as pd
import logging

logging.basicConfig(filename=snakemake.log[0], filemode="w",
                    format="%(asctime)s [%(levelname)s] %(message)s", level=logging.INFO)

logging.info("Merging CEL files by BARCODE")

output_file = snakemake.output.cel
sample_map_file = snakemake.output.sample_map
input_files = snakemake.input.cel

merged_df = None
per_chrom_droplet_ids = {}

for cel_file in input_files:
    chrom = os.path.basename(cel_file).split(".")[0]
    logging.info(f"Processing {chrom} from {cel_file}")

    with gzip.open(cel_file, "rt") as f:
        df = pd.read_csv(f, sep="\t", dtype={
            '#DROPLET_ID': 'Int64',
            'BARCODE': str,
            'NUM.READ': 'Int64',
            'NUM.UMI': 'Int64',
            'NUM.UMIwSNP': 'Int64',
            'NUM.SNP': 'Int64'
        })

    # Record per-chromosome mapping
    per_chrom_droplet_ids[chrom] = df[['BARCODE', '#DROPLET_ID']].set_index('BARCODE')

    # Drop the old ID and prepare for merge
    df = df.drop(columns=['#DROPLET_ID']).set_index('BARCODE')

    if merged_df is None:
        merged_df = df
    else:
        merged_df = merged_df.add(df, fill_value=0)

# Finalize merged data
merged_df = merged_df.reset_index()
merged_df.insert(0, '#DROPLET_ID', range(len(merged_df)))  # new unique IDs

# Save merged CEL file
merged_df.to_csv(output_file, sep='\t', index=False, compression='gzip')
logging.info(f"Saved merged CEL file to {output_file}")

# Build sample map
sample_map = merged_df[['#DROPLET_ID', 'BARCODE']].set_index('BARCODE')
for chrom, mapping in per_chrom_droplet_ids.items():
    sample_map[chrom] = mapping['#DROPLET_ID']

sample_map = sample_map.reset_index()
sample_map.to_csv(sample_map_file, sep='\t', index=False, compression='gzip')
logging.info(f"Saved full sample map to {sample_map_file}")
