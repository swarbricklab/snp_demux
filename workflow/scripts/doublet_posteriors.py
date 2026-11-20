"""Export per-barcode doublet posteriors for one capture.

One row per barcode in the consensus call set (cell_assignment.tsv), joined to
vireo's per-barcode posteriors (donor_ids.tsv). Downstream expression-based
doublet detection (e.g. scDblFinder in mkobj) uses the consensus doublets as
known doublets and the vireo posteriors to estimate how many doublets the hard
calls miss.

Columns:
    capture_id, barcode,
    consensus_status        singlet | recovered | doublet | unassigned (from the consensus in force)
    consensus_assignment    sample_id for singlets, NA otherwise
    vireo_prob_doublet      vireo posterior P(doublet)
    vireo_prob_max          vireo posterior of the best singlet donor
    vireo_best_singlet      vireo's best singlet, as a sample_id of this capture

vireo values are passed through verbatim (no float re-formatting). Barcodes the
consensus has but vireo lacks get NA vireo columns.

Every id column is read as a string: sample ids are often bare digits, and a
reader left to guess types turns "007" into 7.
"""
import sys

import pandas as pd

COLUMNS = [
    "capture_id",
    "barcode",
    "consensus_status",
    "consensus_assignment",
    "vireo_prob_doublet",
    "vireo_prob_max",
    "vireo_best_singlet",
]

# readr writes missing values as "NA"; treat that and empty cells as missing,
# but nothing else (pandas' default list includes e.g. "NaN", "null", "N/A").
NA_VALUES = ["NA", ""]


def read_str_table(path, sep):
    return pd.read_csv(
        path, sep=sep, dtype=str, keep_default_na=False, na_values=NA_VALUES
    )


def donor_to_sample_map(samples, capture_id):
    """donor_id -> sample_id for one capture (first sample wins on duplicates)."""
    cap = samples.loc[samples["capture_id"] == capture_id, ["donor_id", "sample_id"]]
    cap = cap.drop_duplicates(subset="donor_id", keep="first")
    return dict(zip(cap["donor_id"], cap["sample_id"]))


def build_posteriors(cell_assignment, vireo, donor_map, capture_id):
    """Join consensus calls to vireo posteriors; one row per consensus barcode."""
    v = vireo[["cell", "prob_doublet", "prob_max", "best_singlet"]].rename(
        columns={
            "cell": "barcode",
            "prob_doublet": "vireo_prob_doublet",
            "prob_max": "vireo_prob_max",
            "best_singlet": "vireo_best_singlet",
        }
    )
    v = v.drop_duplicates(subset="barcode", keep="first")
    # Sanity check: posteriors must be numeric (fail loudly on a malformed file).
    for col in ("vireo_prob_doublet", "vireo_prob_max"):
        pd.to_numeric(v[col], errors="raise")
    # vireo reports donors; everything under assignment/ is in sample space.
    v["vireo_best_singlet"] = v["vireo_best_singlet"].map(
        lambda d: donor_map.get(d, d) if isinstance(d, str) else d
    )

    c = cell_assignment.rename(
        columns={"status": "consensus_status", "assignment": "consensus_assignment"}
    )[["barcode", "consensus_status", "consensus_assignment"]]

    out = c.merge(v, on="barcode", how="left", validate="one_to_one")
    out.insert(0, "capture_id", capture_id)
    return out[COLUMNS]


def write_csv(df, path):
    df.to_csv(path, index=False, na_rep="NA")


def main(sm):
    sys.stderr = open(sm.log[0], "w")
    capture_id = sm.wildcards.capture_id
    samples = read_str_table(sm.input.samples, ",")
    cell_assignment = read_str_table(sm.input.cell_assignment, "\t")
    vireo = read_str_table(sm.input.vireo, "\t")

    out = build_posteriors(
        cell_assignment, vireo, donor_to_sample_map(samples, capture_id), capture_id
    )
    n_missing = int(out["vireo_prob_doublet"].isna().sum())
    if n_missing:
        print(
            f"{capture_id}: {n_missing} consensus barcodes have no vireo row",
            file=sys.stderr,
        )
    print(out["consensus_status"].value_counts(dropna=False), file=sys.stderr)
    write_csv(out, sm.output.csv)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
