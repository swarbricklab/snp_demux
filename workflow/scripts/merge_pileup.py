import os
import gzip
import pandas as pd
import logging

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    format="%(asctime)s [%(levelname)s] %(message)s",
    level=logging.INFO
)


logging.info("Merging cel files")

# Read individual chromosome cel files into a list of dataframes
cel_merged=None
full_sample_map=None
for cel_file in snakemake.input.cel:
    chromosome=os.path.basename(cel_file).split(".")[0]
    logging.info("Parsing cel file for %s", chromosome)
    with gzip.open(cel_file, "rt") as f:
        cel=pd.read_csv(f,sep='\t').set_index('BARCODE')
    print(cel)
    if cel_merged is None:
        cel_merged = cel  # Initialize the merged DataFrame
    else:
        cel_merged = cel_merged.add(cel, fill_value=0)
    # Extend list of barcode-droplet id maps
    logging.info("Extending sample map with %s", chromosome)
    sample_map=(
        cel.filter(items=['#DROPLET_ID'])
        .rename(columns={'#DROPLET_ID': chromosome})
    )
    if full_sample_map is None:
        full_sample_map = sample_map
    else:
        full_sample_map = pd.merge(
            full_sample_map,
            sample_map,
            how="outer",
            left_index=True,
            right_index=True
        )

logging.info("Wrangling merged cel data")
cel_merged=cel_merged.drop('#DROPLET_ID', axis=1).reset_index()
cel_merged.index.name = '#DROPLET_ID'
logging.info("Saving merged cel file")
cel_merged.to_csv(snakemake.output.cel, sep='\t', index=True, compression='gzip')

# 2. Make make map between the droplet ids used in individual chromosomes and the final (merged) droplet ids
logging.info("Making droplet id map")
full_sample_map=pd.merge(
    cel_merged.reset_index().set_index('BARCODE').filter(items=['#DROPLET_ID']),
    full_sample_map,
    left_index=True,right_index=True
    )
logging.info("Finished droplet id map")
print("Sample Map")
print(full_sample_map)
del(cel_merged)

# 3. merge var files
logging.info("Loading var files")
var_list=[]
for var_file in snakemake.input.var:
    with gzip.open(var_file, "rt") as var:
        df=pd.read_csv(var, sep='\t').rename(columns={'#SNP_ID': 'original_SNP_ID'})
        var_list.append(df)
logging.info("Merging var files for individual chromosomes into a single var file")
var_merged=pd.concat(var_list, ignore_index=True)
var_merged.index.name='#SNP_ID'
(
    var_merged.drop('original_SNP_ID', axis=1)
    .to_csv(snakemake.output.var, sep='\t', index=True, compression='gzip')
)
logging.info("Saved merged var file")
del(var_list)


# 4. Make map of SNP ids between individual chromosomes and merged var file
logging.info("Make SNP map")
snp_map=var_merged.filter(items=['original_SNP_ID','CHROM'])
print("SNP Map")
print(snp_map)
del(var_merged)

# 4. Merge plp files
logging.info("Merging plp files")
plp_list=[]
for plp_file in snakemake.input.plp:
    with gzip.open(plp_file, "rt") as plp:
        chromosome=os.path.basename(plp_file).split(".")[0]
        logging.info("Parsing plp file for %s", chromosome)
        df=(
            pd.read_csv(plp, sep='\t')
            .rename(columns={'SNP_ID': 'original_SNP_ID', '#DROPLET_ID': 'original_DROPLET_ID'})
        )
        print(df)
        snps=(
            snp_map.loc[snp_map['CHROM'] == chromosome]
            .drop('CHROM', axis=1)
            .rename(columns={'#SNP_ID': 'original_SNP_ID'})
            .reset_index()
        )
        print(snps)
        droplets=(
            full_sample_map.filter(items=['#DROPLET_ID', chromosome])
            .rename(columns={chromosome: 'original_DROPLET_ID'})
        )
        print(droplets)
        plp_joined=(
            pd.merge(df,snps, how='left', on='original_SNP_ID')
            .merge(droplets, how='left', on='original_DROPLET_ID')
            .filter(items=['#DROPLET_ID','#SNP_ID','ALLELES','BASEQS'])
            .rename(columns={'#SNP_ID': 'SNP_ID'}) 
        )
        print(plp_joined)
        plp_list.append(plp_joined)

logging.info("Merging plp files for individual chromosomes into a single plp file")
plp_merged=pd.concat(plp_list, ignore_index=True)
print(plp_merged)
del(plp_list)
plp_merged.to_csv(snakemake.output.plp, sep='\t', index=False, compression='gzip')
del(plp_merged)

# 5. Merge umi files
logging.info("Merging umi files")

umi_headers=[
    "original_DROPLET_ID", # Integer ID of the barcode
    "UMI",                 # UMI sequence
    "Total_Fwd_Reads",     # Total length of forward reads
    "Total_Rev_Reads",     # Total length of reverse reads
    "Read_Details"         # Combined read details (optional, will split later)
]
umi_list=[]
for umi_file in snakemake.input.umi:
    chromosome=os.path.basename(umi_file).split(".")[0]
    logging.info("Parsing umi file for %s",chromosome)
    temp_file=umi_file + '_tmp'
    with gzip.open(umi_file, "rt") as umi, open(temp_file, 'wt') as temp:
        for line in umi:
            columns = line.strip().split("\t")
            # Collapse fifth and subsequent columns into a single column
            collapsed = ",".join(columns[4:]) if len(columns) > 4 else ""
            # Write the first four columns followed by the collapsed column
            temp.write("\t".join(columns[:4] + [collapsed]) + "\n")
    with open(temp_file, "rt") as temp:
        df=pd.read_csv(temp, sep='\t', header=None, names=umi_headers)
    os.remove(temp_file)
    droplets=(
        full_sample_map.filter(items=['#DROPLET_ID', chromosome])
        .rename(columns={chromosome: 'original_DROPLET_ID'})
    )
    umi_list.append(
        (
            pd.merge(df,droplets,how='left',on='original_DROPLET_ID')
            .filter(items=['#DROPLET_ID','UMI','Total_Fwd_Reads','Total_Rev_Reads','Read_Details'])
        )
    )

print("Merging umi files for individual chromosomes into a single umi file")
umi_merged=pd.concat(umi_list, ignore_index=True)
#del(umi_list)
umi_merged["Read_Details"] = umi_merged["Read_Details"].str.split(",")
umi_merged=pd.concat(
    [
        umi_merged.drop(columns=["Read_Details"]),
        pd.DataFrame(umi_merged["Read_Details"].to_list(), index=umi_merged.index),
    ], axis=1
)
umi_merged.to_csv(snakemake.output.umi, header=False, index=False, compression='gzip')





