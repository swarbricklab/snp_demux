"""Per-capture doublet summary, and the cohort roll-up of those summaries.

Under random pairing the share of doublets that are cross-donor (and so
visible to genotype demultiplexing) is 1 - sum(p_i^2), where p_i is sample i's
share of the capture's cells. This summary records the observed cross-donor
doublet fraction f_x together with every p_i, so a downstream step can turn
f_x into an expected total doublet rate.

Per-capture columns (one row):
    capture_id, n_barcodes, n_singlet, n_doublet_hard, n_unassigned,
    sum_prob_doublet
        sum of vireo P(doublet) over the capture's barcodes
    doublet_rate_cross_donor
        n_doublet_hard / n_barcodes (f_x, from the hard calls)
    doublet_rate_cross_donor_posterior
        sum_prob_doublet / n_barcodes (f_x, from vireo's posteriors)
    n_samples
        samples declared for the capture in samples.csv
    same_donor_pair_prob
        sum of p_i^2: the chance two randomly paired cells share a donor.
        p_i from singlet counts over ALL declared samples (passengers with
        zero cells included)
    sample_props
        "sample:p;sample:p;..." in samples.csv order

Expected doublet rates (all LOWER BOUNDS: genotype methods miss some
cross-donor doublets, and same-sample aggregates from incomplete dissociation
are always within-donor, so neither is counted here):
    loading_rate_per_1k_cells
        loading-curve doublet rate per 1,000 cells for the capture's chip
        (params.doublet.rate_per_1k)
    doublet_rate_total
        f_x / (1 - same_donor_pair_prob): expected total doublet rate,
        where f_x = doublet_rate_cross_donor
    doublet_rate_total_sd
        binomial SE of f_x, propagated (same_donor_pair_prob treated as fixed):
        sqrt(f_x * (1 - f_x) / n_barcodes) / (1 - same_donor_pair_prob)
    doublet_rate_within_donor
        doublet_rate_total - f_x: expected within-donor (genotype-invisible)
        doublet rate
    doublet_rate_loading
        n_barcodes / 1000 * loading_rate_per_1k_cells: what the chip's loading
        curve predicts
    doublet_rate_vs_loading
        doublet_rate_total / doublet_rate_loading. QC: well above 1 suggests
        overloading, clumping or a demux problem
    doublet_rate_fallback
        why doublet_rate_total, doublet_rate_total_sd and
        doublet_rate_within_donor are NA, else empty:
        single_sample | no_singlets | one_sample_only

Rule `combine_doublet_summary` concatenates the per-capture files.

Full column definitions: docs/data_dictionary.md (keep the two in step).
"""
import math
import sys

import pandas as pd

COLUMNS = [
    "capture_id",
    "n_barcodes",
    "n_singlet",
    "n_doublet_hard",
    "n_unassigned",
    "sum_prob_doublet",
    "doublet_rate_cross_donor",
    "doublet_rate_cross_donor_posterior",
    "n_samples",
    "same_donor_pair_prob",
    "sample_props",
    "loading_rate_per_1k_cells",
    "doublet_rate_total",
    "doublet_rate_total_sd",
    "doublet_rate_within_donor",
    "doublet_rate_loading",
    "doublet_rate_vs_loading",
    "doublet_rate_fallback",
]

NA_VALUES = ["NA", ""]
DIGITS = 6  # rounding for the derived floats written to CSV


def read_str_table(path, sep=","):
    return pd.read_csv(
        path, sep=sep, dtype=str, keep_default_na=False, na_values=NA_VALUES
    )


def declared_samples(samples, capture_id):
    """sample_ids declared for a capture, in samples.csv order."""
    return samples.loc[samples["capture_id"] == capture_id, "sample_id"].tolist()


def sample_proportions(posteriors, samples):
    """p_i for every declared sample, from consensus singlet assignments.

    Returns an ordered dict-like Series (index = declared samples), or None if
    the capture has no singlets assigned to a declared sample.
    """
    assigned = posteriors["consensus_assignment"]
    counts = assigned[assigned.isin(samples)].value_counts()
    counts = counts.reindex(samples, fill_value=0)
    total = counts.sum()
    if total == 0:
        return None
    return counts / total


def _ratio(num, den):
    return num / den if den else float("nan")


def summarise(posteriors, samples, capture_id):
    """One summary row (a dict keyed by COLUMNS) for one capture."""
    status = posteriors["consensus_status"]
    n_barcodes = len(posteriors)
    n_doublet = int((status == "doublet").sum())
    prob_doublet = pd.to_numeric(posteriors["vireo_prob_doublet"], errors="raise")
    sum_prob_doublet = float(prob_doublet.sum(skipna=True))

    props = sample_proportions(posteriors, samples)
    if props is None:
        same_donor_pair_prob = float("nan")
        sample_props = float("nan")
    else:
        same_donor_pair_prob = float((props**2).sum())
        sample_props = ";".join(f"{s}:{round(p, DIGITS)}" for s, p in props.items())

    return {
        "capture_id": capture_id,
        "n_barcodes": n_barcodes,
        "n_singlet": int((status == "singlet").sum()),
        "n_doublet_hard": n_doublet,
        "n_unassigned": int((status == "unassigned").sum()),
        "sum_prob_doublet": sum_prob_doublet,
        "doublet_rate_cross_donor": _ratio(n_doublet, n_barcodes),
        "doublet_rate_cross_donor_posterior": _ratio(sum_prob_doublet, n_barcodes),
        "n_samples": len(samples),
        "same_donor_pair_prob": same_donor_pair_prob,
        "sample_props": sample_props,
    }


def expected_rates(row, rate_per_1k):
    """Expected doublet rates from a summarise() row (see module docstring)."""
    nan = float("nan")
    n = row["n_barcodes"]
    f_x = row["doublet_rate_cross_donor"]
    same_donor_pair_prob = row["same_donor_pair_prob"]
    doublet_rate_loading = n / 1000 * rate_per_1k

    if row["n_samples"] < 2:
        fallback = "single_sample"
    elif math.isnan(same_donor_pair_prob):
        fallback = "no_singlets"
    elif 1 - same_donor_pair_prob <= 0:
        # every singlet came from one sample, so no doublet is expected to be
        # cross-donor and f_x carries no information about doublet_rate_total
        fallback = "one_sample_only"
    else:
        fallback = ""

    if fallback or n == 0:
        doublet_rate_total = doublet_rate_total_sd = doublet_rate_within_donor = nan
    else:
        cross = 1 - same_donor_pair_prob
        doublet_rate_total = f_x / cross
        doublet_rate_total_sd = math.sqrt(f_x * (1 - f_x) / n) / cross
        doublet_rate_within_donor = doublet_rate_total - f_x

    return {
        "loading_rate_per_1k_cells": rate_per_1k,
        "doublet_rate_total": doublet_rate_total,
        "doublet_rate_total_sd": doublet_rate_total_sd,
        "doublet_rate_within_donor": doublet_rate_within_donor,
        "doublet_rate_loading": doublet_rate_loading,
        "doublet_rate_vs_loading": _ratio(doublet_rate_total, doublet_rate_loading),
        "doublet_rate_fallback": fallback,
    }


def to_frame(row, columns=COLUMNS):
    df = pd.DataFrame([row])[columns]
    return df.round(DIGITS)


def write_csv(df, path):
    df.to_csv(path, index=False, na_rep="NA")


def read_calls(cell_assignment, assignments):
    """Per-barcode consensus call + vireo P(doublet), one row per barcode.

    cell_assignment.tsv gives the consensus call (status, sample); the vireo
    posterior comes from sample_barcode_assignments.tsv (Vireo_prob_doublet,
    added by add_posteriors.py).
    """
    calls = cell_assignment.rename(
        columns={"status": "consensus_status", "assignment": "consensus_assignment"}
    )[["barcode", "consensus_status", "consensus_assignment"]]
    vireo = assignments[["Barcode", "Vireo_prob_doublet"]].rename(
        columns={"Barcode": "barcode", "Vireo_prob_doublet": "vireo_prob_doublet"}
    )
    vireo = vireo.drop_duplicates(subset="barcode", keep="first")
    return calls.merge(vireo, on="barcode", how="left", validate="one_to_one")


def main_capture(sm):
    capture_id = sm.wildcards.capture_id
    samples = declared_samples(read_str_table(sm.input.samples), capture_id)
    posteriors = read_calls(
        read_str_table(sm.input.cell_assignment, sep="\t"),
        read_str_table(sm.input.assignments, sep="\t"),
    )
    row = summarise(posteriors, samples, capture_id)
    row.update(expected_rates(row, float(sm.params.rate_per_1k)))
    print(row, file=sys.stderr)
    write_csv(to_frame(row), sm.output.csv)


def main_combine(sm):
    frames = [read_str_table(p) for p in sm.input.summaries]
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    write_csv(combined, sm.output.csv)


def main(sm):
    sys.stderr = open(sm.log[0], "w")
    if sm.rule == "combine_doublet_summary":
        main_combine(sm)
    else:
        main_capture(sm)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
