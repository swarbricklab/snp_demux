#!/usr/bin/env python3

import pandas as pd
import gzip
import argparse
import sys

# ---- CLI Argument Parsing ---- #
parser = argparse.ArgumentParser(description="Sanity check for PLP merge integrity.")
parser.add_argument("--chrom", required=True, help="Chromosome (e.g., chr12)")
parser.add_argument("--droplet", type=int, required=True, help="Original DROPLET_ID (integer)")
parser.add_argument("--snp", type=int, required=True, help="Original SNP_ID (integer)")
parser.add_argument("--capture", required=True, help="Capture ID (e.g., PRADO_9a)")
parser.add_argument("--basedir", default="data/snp/demux/demuxlet/", help="Base directory of input files")
args = parser.parse_args()

chrom = args.chrom
original_droplet_id = args.droplet
original_snp_id = args.snp
capture_id = args.capture
base_dir = args.basedir

# ---- FILE PATHS ---- #
full_sample_map = f"{base_dir}/{capture_id}/full_sample_map.tsv.gz"
snp_map_file = f"{base_dir}/{capture_id}/snp_map.tsv.gz"

plp_orig = f"{base_dir}/pileup/{capture_id}/{chrom}.plp.gz"
cel_orig = f"{base_dir}/pileup/{capture_id}/{chrom}.cel.gz"
var_orig = f"{base_dir}/pileup/{capture_id}/{chrom}.var.gz"
umi_orig = f"{base_dir}/pileup/{capture_id}/{chrom}.umi.gz"

plp_merged = f"{base_dir}/{capture_id}/pileup.plp.gz"
cel_merged = f"{base_dir}/{capture_id}/pileup.cel.gz"
var_merged = f"{base_dir}/{capture_id}/pileup.var.gz"
umi_merged = f"{base_dir}/{capture_id}/pileup.umi.gz"

# ---- LOAD MAPS ---- #
print("Loading sample and SNP maps...")
sample_map = pd.read_csv(full_sample_map, sep="\t", compression="gzip", index_col=0)
snp_map = pd.read_csv(snp_map_file, sep="\t", compression="gzip", index_col=0)

# ---- LOOKUP final IDs ---- #
def get_final_droplet_id(sample_map: pd.DataFrame, chrom: str, original_droplet_id: int) -> int:
    """Return the final merged #DROPLET_ID for a given chromosome and original droplet ID."""
    if chrom not in sample_map.columns:
        raise ValueError(f"Chromosome column '{chrom}' not found in full sample map.")

    matches = sample_map[sample_map[chrom] == original_droplet_id]

    if matches.empty:
        raise ValueError(f"original_DROPLET_ID {original_droplet_id} not found for {chrom} in full sample map.")

    final_id = matches["#DROPLET_ID"].values[0]
    print(f"original_DROPLET_ID {original_droplet_id} → final #DROPLET_ID {final_id}")
    return final_id

def get_final_snp_id(snp_map: pd.DataFrame, chrom: str, original_snp_id: int) -> int:
    """Return the final #SNP_ID (index) for a given chromosome and original SNP ID."""
    matches = snp_map[
        (snp_map["original_SNP_ID"] == original_snp_id) & (snp_map["CHROM"] == chrom)
    ]

    if matches.empty:
        raise ValueError(f"original_SNP_ID {original_snp_id} not found for {chrom} in SNP map.")

    final_id = matches.index.values[0]
    print(f"original_SNP_ID     {original_snp_id} → final #SNP_ID     {final_id}")
    return final_id

# ---- Define helper to search gzipped or plain TSVs ---- #
def grep_tsv(file, col1, val1, col2=None, val2=None):
    with gzip.open(file, "rt") if file.endswith(".gz") else open(file, "rt") as f:
        header = f.readline().strip().split("\t")
        matches = []
        for line in f:
            fields = line.strip().split("\t")
            d = dict(zip(header, fields))
            if d.get(col1) == str(val1) and (col2 is None or d.get(col2) == str(val2)):
                matches.append(d)
        return matches
    
import pandas as pd
import gzip

def lookup_droplet_in_cel(cel_file: str, droplet_id: int) -> pd.Series | None:
    """
    Look up a specific #DROPLET_ID in a .cel.gz or .cel file and return the row as a pandas Series.

    Parameters:
        cel_file (str): Path to .cel.gz or .cel file.
        droplet_id (int): The #DROPLET_ID to search for.

    Returns:
        pd.Series if found, else None.
    """
    open_fn = gzip.open if cel_file.endswith(".gz") else open

    with open_fn(cel_file, "rt") as f:
        df = pd.read_csv(f, sep="\t")

    if "#DROPLET_ID" not in df.columns:
        raise ValueError(f"'#DROPLET_ID' column not found in file: {cel_file}")

    match = df[df["#DROPLET_ID"] == droplet_id]

    if match.empty:
        print(f"#DROPLET_ID {droplet_id} not found in {cel_file}")
        return None

    if len(match) > 1:
        print(f"WARNING: Multiple rows found for #DROPLET_ID {droplet_id} in {cel_file}")

    return match.iloc[0]

def lookup_snp_in_var(var_file: str, snp_id: int) -> pd.Series | None:
    """
    Look up a specific #SNP_ID in a .var.gz or .var file and return the row as a pandas Series.
    """
    import pandas as pd
    import gzip

    open_fn = gzip.open if var_file.endswith(".gz") else open
    with open_fn(var_file, "rt") as f:
        df = pd.read_csv(f, sep="\t", index_col=0)

    if snp_id not in df.index:
        print(f"#SNP_ID {snp_id} not found in {var_file}")
        return None

    return df.loc[snp_id]

def lookup_in_plp(plp_file: str, droplet_id: int, snp_id: int) -> pd.DataFrame:
    """
    Look up all rows for a specific (#DROPLET_ID, SNP_ID) pair in a .plp.gz or .plp file.
    Returns a DataFrame with matching rows.
    """
    import pandas as pd
    import gzip

    open_fn = gzip.open if plp_file.endswith(".gz") else open
    with open_fn(plp_file, "rt") as f:
        df = pd.read_csv(f, sep="\t")
    print(df.shape)

    matches = df[(df["#DROPLET_ID"] == droplet_id) & (df["SNP_ID"] == snp_id)]
    if matches.empty:
        print(f"No match for (#DROPLET_ID={droplet_id}, SNP_ID={snp_id}) in {plp_file}")
    return matches

def lookup_in_umi(umi_file: str, droplet_id: int) -> list[dict]:
    """
    Line-by-line lookup of entries in a non-rectangular .umi(.gz) file for a given original_DROPLET_ID.
    Returns a list of dicts, one per matching row.
    """
    import gzip

    open_fn = gzip.open if umi_file.endswith(".gz") else open
    matches = []

    with open_fn(umi_file, "rt") as f:
        for line in f:
            columns = line.rstrip("\n").split("\t")
            try:
                row_droplet_id = int(columns[0])
            except (IndexError, ValueError):
                continue  # skip malformed lines

            if row_droplet_id == droplet_id:
                matches.append({
                    "original_DROPLET_ID": row_droplet_id,
                    "UMI": columns[1],
                    "Total_Fwd_Reads": int(columns[2]) if len(columns) > 2 else None,
                    "Total_Rev_Reads": int(columns[3]) if len(columns) > 3 else None,
                    "Read_Details": columns[4:] if len(columns) > 4 else []
                })

    if not matches:
        print(f"No UMI matches for original_DROPLET_ID={droplet_id} in {umi_file}")
    return matches


# ----- Main -------

try:
    final_droplet_id = get_final_droplet_id(sample_map, chrom, original_droplet_id)
    final_snp_id = get_final_snp_id(snp_map, chrom, original_snp_id)
except ValueError as e:
    print(f"ERROR: {e}")
    sys.exit(1)




# Compare cel
print(
    lookup_droplet_in_cel(cel_file=cel_orig,droplet_id=original_droplet_id)
)
print(
    lookup_droplet_in_cel(cel_file=cel_merged,droplet_id=final_droplet_id)
)
# Compare var
print(
    lookup_snp_in_var(var_file=var_orig,snp_id=original_snp_id)
)
print(
    lookup_snp_in_var(var_file=var_merged,snp_id=final_snp_id)
)
# Compare plp
print(
    lookup_in_plp(plp_file=plp_orig,droplet_id=original_droplet_id,snp_id=original_snp_id)
)
print(
    lookup_in_plp(plp_file=plp_merged,droplet_id=final_droplet_id,snp_id=final_snp_id)
)
# Compare umi
print(
    lookup_in_umi(umi_file=umi_orig,droplet_id=original_droplet_id)
)
print(
    lookup_in_umi(umi_file=umi_merged,droplet_id=final_droplet_id)
)