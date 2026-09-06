# ===================
# join_sex.py
# -------------------
# Join the per-capture sex_calls.csv (barcode,sex_call,confidence) onto a per-cell
# assignment / barcode table, adding `sex_call` and `sex_confidence` columns.
#
# Robust to a missing or empty sex file: the base table is emitted unchanged plus
# empty sex columns, so downstream (hierarchical consensus) can always read a
# well-formed sex-joined table.
#
# NOTE: snakemake-driven script (reads the `snakemake` object), mirroring
# annotate.py / demuxafy_counts_AtLeastHalfSinglet.py. The sex CSV path arrives via
# params (an optional upstream source with no producing rule), not input.
# ===================

import os
import sys
import logging as log

import pandas as pd

log.basicConfig(
    filename=snakemake.log[0],
    level=log.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s",
)


def _barcode_col(df):
    """Best-effort identification of the barcode column in the base table."""
    for cand in ("barcode", "cell", "cell_barcode"):
        if cand in df.columns:
            return cand
    # souporcell clusters.tsv leads with the barcode column regardless of name.
    return df.columns[0]


base_table = snakemake.input["base_table"]
sex_csv = snakemake.params["sex_csv"]
out_path = snakemake.output["sex_joined"]

log.info("base_table: %s", base_table)
log.info("sex_csv: %s", sex_csv)
log.info("out: %s", out_path)

# souporcell clusters.tsv is tab-separated with a header (barcode, status, ...).
base = pd.read_csv(base_table, sep="\t", dtype=str)
bc_col = _barcode_col(base)
base = base.rename(columns={bc_col: "barcode"})
log.info("base rows: %d, barcode col was: %s", len(base), bc_col)

have_sex = bool(sex_csv) and os.path.exists(sex_csv) and os.path.getsize(sex_csv) > 0

if have_sex:
    try:
        sex = pd.read_csv(sex_csv, dtype=str)
    except Exception as exc:  # empty / malformed -> degrade gracefully
        log.warning("Failed to read sex CSV (%s): %s -- emitting empty sex columns", sex_csv, exc)
        have_sex = False

if have_sex:
    # Tolerate column-name variants for the sex call + its confidence.
    cols = {c.lower(): c for c in sex.columns}
    bc = cols.get("barcode", list(sex.columns)[0])
    call_col = cols.get("sex_call") or cols.get("sex") or cols.get("call")
    conf_col = cols.get("confidence") or cols.get("sex_confidence") or cols.get("conf")

    keep = sex[[bc] + [c for c in (call_col, conf_col) if c]].copy()
    rename = {bc: "barcode"}
    if call_col:
        rename[call_col] = "sex_call"
    if conf_col:
        rename[conf_col] = "sex_confidence"
    keep = keep.rename(columns=rename)
    if "sex_call" not in keep.columns:
        keep["sex_call"] = pd.NA
    if "sex_confidence" not in keep.columns:
        keep["sex_confidence"] = pd.NA
    keep = keep.drop_duplicates(subset="barcode")

    joined = base.merge(keep, on="barcode", how="left")
    log.info("joined %d sex calls onto %d barcodes", keep["sex_call"].notna().sum(), len(base))
else:
    log.info("No usable sex file -- emitting empty sex columns")
    joined = base.copy()
    joined["sex_call"] = pd.NA
    joined["sex_confidence"] = pd.NA

joined.to_csv(out_path, sep="\t", index=False)
log.info("wrote %s", out_path)
