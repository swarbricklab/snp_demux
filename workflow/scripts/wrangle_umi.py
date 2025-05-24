import os
import gzip
import pandas as pd
import logging

logging.basicConfig(filename=snakemake.log[0], filemode="w",
                    format="%(asctime)s [%(levelname)s] %(message)s", level=logging.INFO)

logging.info(f"Processing UMI file for chromosome {snakemake.wildcards.chrom}")

sample_map_file = snakemake.input.sample_map
input_umi_file = snakemake.input.umi
output_umi_file = snakemake.output.umi

umi_headers = ["original_DROPLET_ID", "UMI", "Total_Fwd_Reads", "Total_Rev_Reads", "Read_Details"]

# Load sample map and ensure droplet ID is an integer
logging.info("Loading sample map")
full_sample_map = pd.read_csv(sample_map_file, sep='\t', compression='gzip', index_col=0)
full_sample_map['#DROPLET_ID'] = full_sample_map['#DROPLET_ID'].fillna(0).astype(int)

chromosome = snakemake.wildcards.chrom

with gzip.open(input_umi_file, "rt") as umi, gzip.open(output_umi_file, "wt") as out_f:
    for line in umi:
        columns = line.strip().split("\t")

        # Ensure all rows have 5 columns (merge extra columns into `Read_Details`)
        collapsed = ",".join(columns[4:]) if len(columns) > 4 else ""
        processed_line = "\t".join(columns[:4] + [collapsed])

        # Convert original droplet ID to int
        droplet_id = columns[0]  # original_DROPLET_ID
        try:
            droplet_id = int(droplet_id)
        except ValueError:
            droplet_id = 0  # Handle non-numeric cases

        # Align against full_sample_map
        if chromosome in full_sample_map.columns:
            if droplet_id in full_sample_map[chromosome].values:
                mapped_id = full_sample_map.loc[full_sample_map[chromosome] == droplet_id, '#DROPLET_ID'].values[0]
                processed_line = processed_line.replace(str(droplet_id), str(mapped_id), 1)

        out_f.write(processed_line + "\n")

logging.info(f"Successfully processed UMI file for chromosome {chromosome}")
