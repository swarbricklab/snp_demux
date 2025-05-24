import os
import gzip
import pandas as pd
import logging

logging.basicConfig(filename=snakemake.log[0], filemode="w",
                    format="%(asctime)s [%(levelname)s] %(message)s", level=logging.INFO)

logging.info("Merging umi files")

sample_map_file = snakemake.input.sample_map
output_file = snakemake.output.umi

umi_headers = ["original_DROPLET_ID", "UMI", "Total_Fwd_Reads", "Total_Rev_Reads", "Read_Details"]

# Load sample map
logging.info("Loading sample map")
full_sample_map = pd.read_csv(sample_map_file, sep='\t', compression='gzip', index_col=0)

# **Ensure `original_DROPLET_ID` in `full_sample_map` is int64**
full_sample_map['#DROPLET_ID'] = full_sample_map['#DROPLET_ID'].fillna(0).astype(int)

merged_umi = []  # List to store DataFrames before concatenation

for umi_file in snakemake.input.umi:
    chromosome = os.path.basename(umi_file).split(".")[0]
    logging.info(f"Processing umi file for {chromosome}")

    with gzip.open(umi_file, "rt") as umi:
        lines = [line.strip().split("\t") for line in umi]

    # Process variable-length rows
    for line in lines:
        collapsed = ",".join(line[4:]) if len(line) > 4 else ""
        line[4:] = [collapsed]

    df = pd.DataFrame(lines, columns=umi_headers)

    # **Ensure `original_DROPLET_ID` in `df` is int64**
    df['original_DROPLET_ID'] = pd.to_numeric(df['original_DROPLET_ID'], errors='coerce').fillna(0).astype(int)

    # Merge barcode-droplet map for the current chromosome
    if chromosome in full_sample_map.columns:
        droplets = full_sample_map[['#DROPLET_ID', chromosome]].rename(columns={chromosome: 'original_DROPLET_ID'})
        df = pd.merge(df, droplets, how="left", on="original_DROPLET_ID")

    df["Read_Details"] = df["Read_Details"].str.split(",")

    merged_umi.append(df)  # Store in list to concatenate later

# **Concatenate all processed UMI DataFrames and write final output**
logging.info("Merging all UMI files")
final_umi = pd.concat(merged_umi, ignore_index=True)
final_umi.to_csv(output_file, sep='\t', index=False, compression='gzip')

logging.info("Successfully saved final merged UMI file")
