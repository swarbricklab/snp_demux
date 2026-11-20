#!/usr/bin/env python3
"""Simulate pooled captures with a known doublet rate, as snp_demux outputs.

Each droplet is a singlet (prob 1 - d) or a doublet (prob d) of two cells drawn
independently from the donor proportions p. A genotype demultiplexer sees a
cross-donor doublet as `doublet` (with probability `sensitivity`), but a
within-donor doublet looks like a singlet of that donor. A fraction
`unassigned` of droplets get no call.

So the expected visible (cross-donor) doublet fraction is f_x = d * (1 - sum p^2),
the relation that doublet_summary.py (#60/#61) inverts.

Per capture, the files the doublet rules consume:
  samples.csv                         (shared; one sample per donor, plus optional passengers)
  <capture>/cell_assignment.tsv        barcode, assignment, status  (current snp_demux format)
  <capture>/cell_assignment.legacy.tsv status NA for every non-singlet (pre-#22 format,
                                       as seen in older mkobj demo data)
  <capture>/donor_ids.tsv              vireo-like posteriors (prob_doublet high for doublets)
  truth.tsv                            the true parameters per capture

Usage: make_doublet_captures.py OUTDIR [--seed 1]
The scenario grid is defined in SCENARIOS below.
"""
import argparse
import csv
import random
from pathlib import Path

# name: (n_droplets, true d, donor proportions, sensitivity, unassigned frac, n_passengers)
SCENARIOS = {
    "equal4_6pct": (5000, 0.08, [0.25] * 4, 1.0, 0.0, 0),  # #61 acceptance case: f_x = 0.06
    "equal8": (8000, 0.10, [1 / 8] * 8, 1.0, 0.01, 0),
    "skewed4": (6000, 0.08, [0.55, 0.25, 0.15, 0.05], 1.0, 0.01, 0),
    "passenger": (4000, 0.06, [0.5, 0.3, 0.2], 1.0, 0.0, 2),  # 2 declared samples, 0 cells
    "two_donor": (4000, 0.08, [0.5, 0.5], 1.0, 0.0, 0),
    "sens80": (6000, 0.08, [0.25] * 4, 0.8, 0.0, 0),  # misses 20% of cross-donor doublets
    "single": (3000, 0.05, [1.0], 1.0, 0.0, 0),  # single-sample capture -> fallback
}


def bc(i):
    return f"BC{i:07d}-1"


def simulate(rng, n, d, p, sens, unassigned, donors):
    rows = []
    for i in range(n):
        if rng.random() < unassigned:
            rows.append((bc(i), "unassigned", None, rng.choice(donors), rng.uniform(0, 0.5)))
            continue
        a = rng.choices(donors, weights=p)[0]
        if rng.random() < d:
            b = rng.choices(donors, weights=p)[0]
            if a != b and rng.random() < sens:
                rows.append((bc(i), "doublet", None, a, rng.uniform(0.6, 1.0)))
            else:  # within-donor (or missed) doublet: indistinguishable from a singlet
                rows.append((bc(i), "singlet", a, a, rng.uniform(0.0, 0.4)))
        else:
            rows.append((bc(i), "singlet", a, a, rng.uniform(0.0, 0.05)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    a.outdir.mkdir(parents=True, exist_ok=True)
    samples_rows, truth = [], []
    for cap, (n, d, p, sens, unas, n_pass) in SCENARIOS.items():
        donors = [f"{cap}_{k:03d}" for k in range(len(p))]
        samples_rows += [(dn, dn, cap) for dn in donors]
        samples_rows += [(f"{cap}_p{k}", f"{cap}_p{k}", cap) for k in range(n_pass)]
        rows = simulate(rng, n, d, p, sens, unas, donors)
        cdir = a.outdir / cap
        cdir.mkdir(exist_ok=True)
        with open(cdir / "cell_assignment.tsv", "w") as f:
            f.write("barcode\tassignment\tstatus\n")
            for b, st, asg, _, _ in rows:
                f.write(f"{b}\t{asg or 'NA'}\t{st}\n")
        with open(cdir / "cell_assignment.legacy.tsv", "w") as f:
            f.write("barcode\tassignment\tstatus\n")
            for b, st, asg, _, _ in rows:
                f.write(f"{b}\t{asg or 'NA'}\t{st if st == 'singlet' else 'NA'}\n")
        with open(cdir / "donor_ids.tsv", "w") as f:
            f.write("cell\tdonor_id\tprob_max\tprob_doublet\tn_vars\tbest_singlet\tbest_doublet\tdoublet_logLikRatio\n")
            for b, st, asg, best, pd_ in rows:
                f.write(f"{b}\t{asg or st}\t{1 - pd_:.2e}\t{pd_:.2e}\t100\t{best}\t{best},{best}\t0.0\n")
        sum_p_sq = sum(x * x for x in p)
        truth.append(dict(capture_id=cap, n=n, d_true=d, sum_p_sq_true=round(sum_p_sq, 6),
                          f_x_expected=round(d * (1 - sum_p_sq) * sens * (1 - unas), 6),
                          sensitivity=sens, unassigned=unas, n_samples=len(p) + n_pass))
    with open(a.outdir / "samples.csv", "w") as f:
        f.write("sample_id,donor_id,capture_id\n")
        f.writelines(f"{s},{dn},{c}\n" for s, dn, c in samples_rows)
    with open(a.outdir / "truth.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(truth[0]), delimiter="\t")
        w.writeheader()
        w.writerows(truth)


if __name__ == "__main__":
    main()
