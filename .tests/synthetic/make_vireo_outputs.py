#!/usr/bin/env python3
"""Generate synthetic vireo outputs for exercising relabel_vireo.R.

Writes one directory per scenario, each holding the files the `relabel_vireo`
rule consumes (donor_ids.tsv, prob_doublet.tsv.gz, prob_singlet.tsv.gz,
prop_ambient.tsv, summary.tsv, two PDFs) plus a samples.csv for one capture.

Scenarios
  missing1     one pool donor has no genotype -> vireo labels it `donor0`
  all_known    every donor genotyped (relabel_vireo.R early-exit path)
  zero_cells   missing1, plus a genotyped donor that receives no cells
               (it appears in the prob_* headers but not in donor_ids.tsv)
  missing2     two donors without genotype (`donor0`, `donor1`) -> script must stop()

Formats mirror real vireo 0.5.x output (see .tests/integration/golden/vireo/).
The PDF figures are written by the driver script in R (a real one-page PDF for
fig_GT_distance_estimated and a 0-byte placeholder for fig_GT_distance_input),
so this generator needs only the Python standard library.

Usage: make_vireo_outputs.py OUTDIR [--cells 25000] [--seed 1]
"""
import argparse
import gzip
import itertools
import random
from pathlib import Path

CAPTURE = "synth_cap"
DONORS = ["101_101", "202_202", "303_303", "404_404"]  # bare-digit-style ids, like OneK1K


def barcode(rng):
    return "".join(rng.choice("ACGT") for _ in range(16)) + "-1"


def scenario(outdir: Path, name: str, labels: list, assigned: list, n_cells: int, rng):
    """labels: donor labels as vireo emits them (genotyped ids or donorN);
    assigned: the subset of labels that actually receive cells."""
    d = outdir / name
    d.mkdir(parents=True, exist_ok=True)
    # samples.csv: the *real* donor ids for the capture (one row per sample)
    with open(d / "samples.csv", "w") as f:
        f.write("sample_id,donor_id,capture_id\n")
        for donor in DONORS:
            f.write(f"{donor},{donor},{CAPTURE}\n")

    pairs = list(itertools.combinations(labels, 2))
    seen = set()
    cells = []
    while len(cells) < n_cells:
        bc = barcode(rng)
        if bc not in seen:
            seen.add(bc)
            cells.append(bc)

    counts = {}
    with open(d / "donor_ids.tsv", "w") as f:
        f.write("cell\tdonor_id\tprob_max\tprob_doublet\tn_vars\tbest_singlet\tbest_doublet\tdoublet_logLikRatio\n")
        for bc in cells:
            r = rng.random()
            best = rng.choice(assigned)
            dbl = rng.choice([p for p in pairs if p[0] in assigned and p[1] in assigned])
            call = best if r < 0.85 else ("doublet" if r < 0.95 else "unassigned")
            counts[call] = counts.get(call, 0) + 1
            f.write(f"{bc}\t{call}\t{rng.random():.2e}\t{rng.random():.2e}\t{rng.randint(5, 400)}"
                    f"\t{best}\t{dbl[0]},{dbl[1]}\t{rng.uniform(-20, 5):.3f}\n")

    with gzip.open(d / "prob_doublet.tsv.gz", "wt") as f:
        f.write("cell\t" + "\t".join(f"{a},{b}" for a, b in pairs) + "\n")
        for bc in cells:
            f.write(bc + "\t" + "\t".join(f"{rng.random():.2e}" for _ in pairs) + "\n")
    with gzip.open(d / "prob_singlet.tsv.gz", "wt") as f:
        f.write("cell\t" + "\t".join(labels) + "\n")
        for bc in cells:
            f.write(bc + "\t" + "\t".join(f"{rng.random():.2e}" for _ in labels) + "\n")
    with open(d / "prop_ambient.tsv", "w") as f:
        f.write("cell\t" + "\t".join(labels) + "\tlogLik_ratio\n")
        for bc in cells:
            f.write(bc + "\t" + "\t".join(f"{rng.random():.4e}" for _ in labels) + f"\t{rng.uniform(-3, 3):.2f}\n")
    with open(d / "summary.tsv", "w") as f:
        f.write("Var1\tFreq\n")
        for k in sorted(counts):
            f.write(f"{k}\t{counts[k]}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--cells", type=int, default=25000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    # vireo sorts donor columns: genotyped ids first, then donor0.. for the learned ones
    m1 = DONORS[:3] + ["donor0"]
    scenario(a.outdir, "missing1", m1, m1, a.cells, rng)
    scenario(a.outdir, "all_known", DONORS, DONORS, a.cells, rng)
    scenario(a.outdir, "zero_cells", m1, [x for x in m1 if x != "202_202"], a.cells, rng)
    m2 = DONORS[:2] + ["donor0", "donor1"]
    scenario(a.outdir, "missing2", m2, m2, a.cells, rng)


if __name__ == "__main__":
    main()
