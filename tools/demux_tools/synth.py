"""Synthetic captures for tests and experiments; the `synth` command.

A synthetic capture is a cellsnp-lite pileup drawn from known donor genotypes,
with a donor VCF and a truth table. Reads follow vireo's model: at a site where
a cell's donor has genotype g (0, 1 or 2 alternative alleles), each UMI is
alternative with probability theta[g]. A doublet of donors a and b with share s
draws each UMI from a with probability s, so its alternative probability is
s * theta[g_a] + (1 - s) * theta[g_b].

Random numbers come from numpy's legacy RandomState, whose stream is fixed
across numpy versions, so a seed gives the same capture in every supported
environment. The draws are made in a fixed order, documented in `make_capture`.
"""
import os
from collections import namedtuple

import numpy as np
import pandas as pd

from demux_tools._io import UsageError, open_text, start_clock, write_provenance, write_table
from demux_tools.cellsnp import Pileup, write_cellsnp

BASES = np.array(list("ACGT"), dtype=object)
GT_TEXT = {0: "0/0", 1: "0/1", 2: "1/1"}

SynthCapture = namedtuple("SynthCapture", ["donors", "pileup", "genotypes", "blank", "truth", "theta", "seed"])
SynthCapture.__doc__ = """A synthetic capture in memory.

donors:    list of donor names, in genotype-column order.
pileup:    Pileup with layers AD, DP and OTH (and its chrom, pos, ref, alt, barcodes).
genotypes: int8 array (n_sites, n_donors), the true alternative-allele count 0, 1 or 2.
blank:     bool array (n_sites, n_donors), genotypes written as missing in the VCF
           (the reads still follow the true genotype).
truth:     pandas DataFrame, one row per barcode in pileup order, all text:
           barcode, kind ("singlet" | "doublet"), donor (singlets; "" for
           doublets), pair (doublets: pair_key of the two donors; "" for
           singlets), donor_a (the singlet's donor, or a doublet's first
           donor), share_a (share of UMIs from donor_a, "%g").
theta:     the three alternative-allele probabilities used.
seed:      the seed.
"""


def default_donor_names(n_donors):
    """"d0", "d1", ... zero-padded so string order equals index order (d00..d11 for 12)."""
    width = len(str(max(n_donors - 1, 0)))
    return ["d%0*d" % (width, k) for k in range(n_donors)]


def make_capture(n_donors=3, n_sites=20, n_cells=10, seed=0, theta=(0.01, 0.5, 0.99), depth=(1, 4),
                 coverage=0.6, doublets=((0, 1, 0.5),), blank_rate=0.0, donor_names=None, oth_rate=0.0,
                 n_chroms=2):
    """Draw a synthetic capture.

    n_donors, n_sites: genotyped donors and piled sites.
    n_cells: singlets; singlet c belongs to donor c % n_donors.
    doublets: one extra cell per (a, b, share) entry, appended after the
        singlets: donors a and b by index, `share` of UMIs from a.
    theta: alternative-allele probability for genotypes 0, 1, 2.
    depth: (low, high) UMIs at a covered site, inclusive.
    coverage: probability that a cell has UMIs at a site.
    blank_rate: probability that a donor's genotype at a site is written as
        missing ("./.") in the VCF.
    donor_names: names in index order (default `default_donor_names`).
    oth_rate: probability that a covered entry also has one UMI of another
        allele (OTH layer).
    n_chroms: sites are spread over chr1..chrN in blocks, positions increasing.

    Every site is polymorphic among the donors when n_donors >= 2. Draw order:
    positions, REF, ALT shift, genotypes, blanks, barcodes, then per entry
    coverage, depth, alternative UMIs and OTH. Returns a SynthCapture.
    """
    if n_donors < 1 or n_sites < 1 or n_cells < 0:
        raise ValueError("need at least one donor and one site")
    lo, hi = int(depth[0]), int(depth[1])
    if lo < 1 or hi < lo:
        raise ValueError("depth must be (low, high) with 1 <= low <= high")
    donors = list(donor_names) if donor_names is not None else default_donor_names(n_donors)
    if len(donors) != n_donors or len(set(donors)) != n_donors:
        raise ValueError("donor_names must be %d distinct names" % n_donors)
    doublets = [(int(a), int(b), float(s)) for a, b, s in doublets]
    for a, b, s in doublets:
        if not (0 <= a < n_donors and 0 <= b < n_donors and a != b and 0.0 <= s <= 1.0):
            raise ValueError("bad doublet (%d, %d, %g)" % (a, b, s))
    theta = tuple(float(t) for t in theta)
    rng = np.random.RandomState(seed)

    # Sites
    n_chroms = max(1, min(int(n_chroms), n_sites))
    chrom_index = (np.arange(n_sites) * n_chroms) // n_sites
    within = np.arange(n_sites) - np.searchsorted(chrom_index, chrom_index)
    pos = 10000 + 1000 * within + rng.randint(0, 500, size=n_sites)
    chrom = np.array(["chr%d" % (c + 1) for c in chrom_index], dtype=object)
    ref_code = rng.randint(0, 4, size=n_sites)
    alt_code = (ref_code + 1 + rng.randint(0, 3, size=n_sites)) % 4

    # Genotypes and blanks
    gt = rng.randint(0, 3, size=(n_sites, n_donors)).astype(np.int8)
    if n_donors >= 2:
        for s in np.flatnonzero((gt == gt[:, :1]).all(axis=1)):
            d = s % n_donors
            gt[s, d] = (gt[s, d] + 1) % 3
    blank = rng.random_sample((n_sites, n_donors)) < blank_rate

    # Cells
    n_total = n_cells + len(doublets)
    barcodes = []
    seen = set()
    while len(barcodes) < n_total:
        bc = "".join(BASES[rng.randint(0, 4, size=16)]) + "-1"
        if bc not in seen:
            seen.add(bc)
            barcodes.append(bc)
    p_alt = np.zeros((n_sites, n_total), dtype=np.float64)
    th = np.asarray(theta)
    rows = []
    for c in range(n_cells):
        d = c % n_donors
        p_alt[:, c] = th[gt[:, d]]
        rows.append((barcodes[c], "singlet", donors[d], "", donors[d], "1"))
    for k, (a, b, s) in enumerate(doublets):
        c = n_cells + k
        p_alt[:, c] = s * th[gt[:, a]] + (1.0 - s) * th[gt[:, b]]
        pair = ",".join(sorted((donors[a], donors[b])))
        rows.append((barcodes[c], "doublet", "", pair, donors[a], "%g" % s))

    # Reads
    covered = rng.random_sample((n_sites, n_total)) < coverage
    dp = rng.randint(lo, hi + 1, size=(n_sites, n_total)) * covered
    ad = rng.binomial(dp, p_alt)
    oth = (rng.random_sample((n_sites, n_total)) < oth_rate) & covered

    def coo(m):
        r, c = np.nonzero(m)
        return r, c, m[r, c]

    pileup = Pileup(chrom, pos, BASES[ref_code], BASES[alt_code], barcodes,
                    layers={"AD": coo(ad), "DP": coo(dp), "OTH": coo(oth.astype(np.int64))})
    truth = pd.DataFrame(rows, columns=["barcode", "kind", "donor", "pair", "donor_a", "share_a"])
    return SynthCapture(donors, pileup, gt, blank, truth, theta, seed)


def _gp_values(g, gp_noise):
    p = np.full(3, gp_noise / 3.0)
    p[g] += 1.0 - gp_noise
    return p


def _pl_text(p):
    with np.errstate(divide="ignore"):
        pl = np.minimum(255, np.round(-10.0 * np.log10(p))).astype(int)
    return ",".join(str(x) for x in pl - pl.min())


def write_vcf(sc, path, tags=("GT",), gp_noise=0.0, rows=None):
    """Write the donor genotypes of `sc` as a VCF.

    tags: FORMAT fields, any of "GT", "GP", "PL", in the order given.
        GT: "0/0", "0/1", "1/1"; "./." when blank.
        GP: (1 - gp_noise) on the true genotype plus gp_noise / 3 on each, "%g"; "." when blank.
        PL: -10 log10(GP) rounded, capped at 255, minus its minimum; "." when blank.
    rows: 0-based site indices to write (default all), e.g. to leave sites out
        of the panel; written in position order as given.
    INFO holds AF, the alternative-allele frequency over non-blank donors ("%.4g").
    Plain gzip with mtime 0 for ".gz" (not BGZF). Returns `path`.
    """
    tags = list(tags)
    bad = [t for t in tags if t not in ("GT", "GP", "PL")]
    if bad or not tags:
        raise ValueError("tags must be from GT, GP, PL; got %s" % tags)
    p = sc.pileup
    rows = np.arange(p.n_variants) if rows is None else np.asarray(rows, dtype=np.int64)
    header = ["##fileformat=VCFv4.2"]
    for c in dict.fromkeys(p.chrom[rows]):
        header.append("##contig=<ID=%s>" % c)
    header.append('##INFO=<ID=AF,Number=A,Type=Float,Description="Alternative allele frequency among the donors">')
    descriptions = {
        "GT": '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "GP": '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype probabilities">',
        "PL": '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">',
    }
    header.extend(descriptions[t] for t in tags)
    header.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + list(sc.donors)))
    with open_text(path, "wt") as f:
        f.write("\n".join(header) + "\n")
        for s in rows:
            g_row, b_row = sc.genotypes[s], sc.blank[s]
            known = g_row[~b_row]
            af = "%.4g" % (known.sum() / (2.0 * len(known))) if len(known) else "."
            fields = [p.chrom[s], str(p.pos[s]), ".", p.ref[s], p.alt[s], ".", "PASS",
                      "AF=%s" % af, ":".join(tags)]
            for g, b in zip(g_row, b_row):
                values = []
                for t in tags:
                    if b:
                        values.append("./." if t == "GT" else ".")
                    elif t == "GT":
                        values.append(GT_TEXT[int(g)])
                    elif t == "GP":
                        values.append(",".join("%g" % x for x in _gp_values(int(g), gp_noise)))
                    else:
                        values.append(_pl_text(_gp_values(int(g), gp_noise)))
                fields.append(":".join(values))
            f.write("\t".join(fields) + "\n")
    return os.fspath(path)


def write_capture(sc, out_dir, tags=("GT",), gp_noise=0.0, vcf_rows=None):
    """Write a synthetic capture to `out_dir` (created if needed).

    Writes out_dir/cellsnp/ (cellsnp-lite files, layers AD, DP, OTH),
    out_dir/donors.vcf.gz (`write_vcf` with tags, gp_noise, rows=vcf_rows) and
    out_dir/truth.tsv. Returns {"dir", "cellsnp", "vcf", "truth", "donors"}.
    The same capture gives the same bytes on every call.
    """
    out_dir = os.fspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    cellsnp_dir = os.path.join(out_dir, "cellsnp")
    write_cellsnp(sc.pileup, cellsnp_dir)
    vcf = write_vcf(sc, os.path.join(out_dir, "donors.vcf.gz"), tags=tags, gp_noise=gp_noise, rows=vcf_rows)
    truth = write_table(sc.truth, os.path.join(out_dir, "truth.tsv"))
    return {"dir": out_dir, "cellsnp": cellsnp_dir, "vcf": vcf, "truth": truth, "donors": list(sc.donors)}


def build_on_disk(out_dir, tags=("GT",), gp_noise=0.0, vcf_rows=None, **make_kwargs):
    """`make_capture(**make_kwargs)` then `write_capture`; the returned dict also holds "capture"."""
    sc = make_capture(**make_kwargs)
    paths = write_capture(sc, out_dir, tags=tags, gp_noise=gp_noise, vcf_rows=vcf_rows)
    paths["capture"] = sc
    return paths


# ---------------------------------------------------------------------------
# Command: synth


def _pair_range(text, cast, what):
    parts = text.split(":")
    if len(parts) != 2:
        raise UsageError("%s must be LOW:HIGH, not %r" % (what, text))
    return cast(parts[0]), cast(parts[1])


def _doublet(text):
    parts = text.split(":")
    if len(parts) != 3:
        raise UsageError("--doublet must be A:B:SHARE (donor indices and A's share), not %r" % text)
    return int(parts[0]), int(parts[1]), float(parts[2])


def configure_synth(parser):
    parser.description = (
        "Write a synthetic capture: OUT/cellsnp/ (cellsnp-lite files), OUT/donors.vcf.gz, "
        "OUT/truth.tsv and OUT/provenance.json. The same options give the same files.")
    parser.add_argument("--out", required=True, help="output folder (created if needed)")
    parser.add_argument("--donors", type=int, default=3, help="genotyped donors (default 3)")
    parser.add_argument("--sites", type=int, default=20, help="piled sites (default 20)")
    parser.add_argument("--cells", type=int, default=10, help="singlet cells (default 10)")
    parser.add_argument("--seed", type=int, default=0, help="random seed (default 0)")
    parser.add_argument("--doublet", action="append", default=None, metavar="A:B:SHARE",
                        help="add a doublet of donors A and B (0-based) with A's share of UMIs; "
                             "repeatable (default 0:1:0.5 when there are two or more donors)")
    parser.add_argument("--no-doublets", action="store_true", help="singlets only")
    parser.add_argument("--blank-rate", type=float, default=0.0, help="share of genotypes written ./. (default 0)")
    parser.add_argument("--coverage", type=float, default=0.6, help="chance a cell has UMIs at a site (default 0.6)")
    parser.add_argument("--depth", default="1:4", help="LOW:HIGH UMIs at a covered site (default 1:4)")
    parser.add_argument("--theta", default="0.01,0.5,0.99", help="alt probability for 0/0,0/1,1/1")
    parser.add_argument("--oth-rate", type=float, default=0.0, help="chance of one OTH UMI per covered entry")
    parser.add_argument("--tags", default="GT", help="VCF FORMAT fields, from GT,GP,PL (default GT)")
    parser.add_argument("--gp-noise", type=float, default=0.0, help="probability mass spread over genotypes in GP/PL")
    parser.add_argument("--donor-names", default=None, help="comma-separated donor names (default d0,d1,...)")


def run_synth(args):
    clock = start_clock()
    if args.no_doublets:
        doublets = []
    elif args.doublet:
        doublets = [_doublet(t) for t in args.doublet]
    else:
        doublets = [(0, 1, 0.5)] if args.donors >= 2 else []
    theta = tuple(float(x) for x in args.theta.split(","))
    if len(theta) != 3:
        raise UsageError("--theta needs three values")
    names = args.donor_names.split(",") if args.donor_names else None
    try:
        sc = make_capture(n_donors=args.donors, n_sites=args.sites, n_cells=args.cells, seed=args.seed,
                          theta=theta, depth=_pair_range(args.depth, int, "--depth"), coverage=args.coverage,
                          doublets=doublets, blank_rate=args.blank_rate, donor_names=names,
                          oth_rate=args.oth_rate)
        paths = write_capture(sc, args.out, tags=args.tags.split(","), gp_noise=args.gp_noise)
    except ValueError as e:
        raise UsageError(str(e))
    summary = {"donors": sc.donors, "sites": sc.pileup.n_variants, "cells": sc.pileup.n_barcodes,
               "doublets": len(doublets), "entries": {k: int(len(v[0])) for k, v in sc.pileup.layers.items()}}
    write_provenance(args.out, args.argv, [], clock=clock, extra=summary)
    print("synth: %d donors, %d sites, %d cells (%d doublets) in %s"
          % (len(sc.donors), sc.pileup.n_variants, sc.pileup.n_barcodes, len(doublets), paths["dir"]))
    return 0
