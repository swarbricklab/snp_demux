#!/usr/bin/env python3
"""Differential test: workflow/scripts/consensus.py vs Demuxafy's Combine_Results.R.

For each seed, simulate one capture's four demultiplexer outputs plus the vireoSNP
souporcell annotation, then build the final assignment tables two ways:

  old  Combine_Results.R --method AtLeastHalfSinglet --ref Vireo   (Demuxafy image)
       -> add_vireosnp_to_combined_results.R (from OLD_REF)          (R image)
       -> wrangle_assignments.R                                      (R image)
  new  consensus.py -> wrangle_assignments.R

and compare sample_barcode_assignments.tsv and cell_assignment.tsv byte for byte.

The simulation is built to reach the corners of Demuxafy's logic: few barcodes
(so agreement proportions tie, including at pct_agreement <= 0.5 where one
cluster gets several candidate donors), clusters that fail to map, junk
clusters, barcodes missing from or extra to each method, mostly-ambiguous
demuxlet, and donor-identity ties between methods.

Two Demuxafy failures are classified rather than compared (the port does not
reproduce them): Combine_Results.R crashes when no cluster of a method reaches
pct_agreement, or when duplicate cluster mappings (pct_agreement <= 0.5) outgrow
data.table's cartesian-join limit ("demuxafy-crash"), and it loses every consensus identity when each
barcode has exactly one distinct donor name, via R's apply() simplification
("demuxafy-bug"). Anything else that differs is a porting error.

Usage: consensus_vs_demuxafy.py [--seeds 1-60] [--jobs 4] [--old-ref origin/main] [--workdir DIR]
Needs singularity and the cached snp_demux-demuxafy and snp_demux-r images.
"""
import argparse
import concurrent.futures as cf
import gzip
import hashlib
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
IMG_DIR = Path("/g/data/a56/software/singularity/images")


def image(url):
    return str(IMG_DIR / (hashlib.md5(url.encode()).hexdigest() + ".simg"))


DEMUXAFY_IMG = image("docker://ghcr.io/swarbricklab/snp_demux-demuxafy:63fb228")
R_IMG = image("docker://ghcr.io/swarbricklab/snp_demux-r:r4.2")
POPSCLE_HEADER = "INT_ID\tBARCODE\tNUM.SNPS\tNUM.READS\tDROPLET.TYPE\tBEST.GUESS\n"


def bc(rng):
    return "".join(rng.choice("ACGT") for _ in range(12)) + "-1"


def simulate(seed, d):
    rng = random.Random(seed)
    n_donors = rng.randint(2, 6)
    donors = [rng.choice(["{0}_{0}", "D{0}"]).format(100 + i) for i in range(n_donors)]
    n = rng.choice([12, 30, 60, 150, 400])
    pct = rng.choice([0.5, 0.7, 0.7, round(rng.uniform(0.3, 0.9), 2)])
    noise = rng.choice([0.0, 0.05, 0.2, 0.5])
    p_amb = {m: rng.choice([0.0, 0.05, 0.3, 0.8]) for m in ("dem", "fmx", "sop", "vir")}
    barcodes = []
    while len(barcodes) < n:
        b = bc(rng)
        if b not in barcodes:
            barcodes.append(b)
    truth = {b: (rng.choice(donors) if rng.random() > 0.1 else "DBL") for b in barcodes}
    junk = rng.random() < 0.3  # an extra cluster soaking up random cells

    def subset(m):  # each method loses some barcodes and gains a few extras
        keep = [b for b in barcodes if rng.random() > rng.choice([0, 0, 0.05, 0.2])]
        extra = [bc(rng) for _ in range(rng.choice([0, 0, 2]))]
        out = keep + extra
        if m != "dem":
            rng.shuffle(out)
        return out

    def call(b, m, labels):
        """(type, first, second) for barcode b, labels = donor -> label."""
        t = truth.get(b, rng.choice(donors))
        if rng.random() < p_amb[m]:
            return "AMB", rng.choice(list(labels.values())), rng.choice(list(labels.values()))
        if t == "DBL" or rng.random() < 0.03:
            a, c = rng.sample(list(labels.values()), 2) if len(labels) > 1 else (labels[donors[0]],) * 2
            return "DBL", a, c
        lab = labels[t] if rng.random() > noise else rng.choice(list(labels.values()))
        return "SNG", lab, lab

    def clusters():
        ids = [str(i) for i in range(n_donors + (1 if junk else 0))]
        rng.shuffle(ids)
        lab = dict(zip(donors, ids))
        return lab, (ids[-1] if junk else None)

    for sub in ("demuxlet", "freemuxlet", "souporcell", "vireo", "vireosnp"):
        (d / sub).mkdir(parents=True, exist_ok=True)

    ident = {x: x for x in donors}
    with open(d / "demuxlet/demuxlet.best", "w") as f:
        f.write(POPSCLE_HEADER)
        for i, b in enumerate(subset("dem")):
            t, a, c = call(b, "dem", ident)
            f.write("%d\t%s\t10\t10\t%s\t%s,%s,%s\n" % (i, b, t, a, c, "0.50" if t == "DBL" else "0.00"))
    lab, junk_id = clusters()
    with gzip.open(d / "freemuxlet/freemuxlet.clust1.samples.gz", "wt") as f:
        f.write(POPSCLE_HEADER)
        for i, b in enumerate(subset("fmx")):
            t, a, c = call(b, "fmx", lab)
            if junk_id and t == "SNG" and rng.random() < 0.1:
                a = c = junk_id
            f.write("%d\t%s\t10\t10\t%s\t%s,%s\n" % (i, b, t, a, c))
    lab, junk_id = clusters()
    with open(d / "souporcell/clusters.tsv", "w") as f, open(d / "vireosnp/clusters_annotated.tsv", "w") as g:
        f.write("barcode\tstatus\tassignment\n")
        g.write("barcode\tstatus\tassignment\n")
        inv = {v: k for k, v in lab.items()}
        for b in subset("sop"):
            t, a, c = call(b, "sop", lab)
            if junk_id and t == "SNG" and rng.random() < 0.1:
                a = c = junk_id
            status = {"SNG": "singlet", "DBL": "doublet", "AMB": "unassigned"}[t]
            asg = a if t != "DBL" else "%s/%s" % (a, c)
            f.write("%s\t%s\t%s\n" % (b, status, asg))
            g.write("%s\t%s\t%s\n" % (b, status, inv.get(asg, asg) if t == "SNG" else asg))
    with open(d / "vireo/donor_ids.tsv", "w") as f:
        f.write("cell\tdonor_id\tprob_max\n")
        for b in subset("vir"):
            t, a, _ = call(b, "vir", ident)
            f.write("%s\t%s\t0.9\n" % (b, {"SNG": a, "DBL": "doublet", "AMB": "unassigned"}[t]))
    with open(d / "samples.csv", "w") as f:
        f.write("sample_id,donor_id,capture_id\n")
        for x in donors:
            f.write("S%s,%s,cap\n" % (x, x))
    return pct


MOCK = """setClass("Snakemake", representation(input = "list", output = "list", params = "list", wildcards = "list", log = "list"))
snakemake <- new("Snakemake", input = list({inp}), output = list({out}), params = list({par}),
  wildcards = list(capture_id = "cap"), log = list("{log}"))
source("{script}")
"""


def run(cmd, tries=3, **kw):
    """subprocess.run with retries: singularity occasionally fails to start
    (exit 255) when several containers launch at once."""
    for i in range(tries):
        try:
            return subprocess.run(cmd, check=True, capture_output=True, **kw)
        except subprocess.CalledProcessError:
            if i == tries - 1:
                raise


def rscript(d, name, script, inp, out, par=""):
    (d / (name + ".R")).write_text(MOCK.format(inp=inp, out=out, par=par, log=d / (name + ".log"), script=script))
    run(["singularity", "exec", "--bind", "%s,%s" % (REPO, d), R_IMG, "Rscript", str(d / (name + ".R"))])


def wrangle(d, calls, outdir):
    outdir.mkdir(exist_ok=True)
    rscript(d, "wrangle_" + outdir.name, REPO / "workflow/scripts/wrangle_assignments.R",
            'sample_map = "%s", demuxafy_assignments = "%s"' % (d / "samples.csv", calls),
            'sample_barcode_assignments = "%s", cell_assignment = "%s"'
            % (outdir / "sample_barcode_assignments.tsv", outdir / "cell_assignment.tsv"),
            'assignment_column = "AtLeastHalfSinglet_Individual_Assignment", '
            'droplettype_column = "AtLeastHalfSinglet_DropletType"')


def run_seed(seed, root, old_add_vireosnp):
    d = root / ("s%d" % seed)
    shutil.rmtree(d, ignore_errors=True)
    pct = simulate(seed, d)
    # old: Demuxafy
    try:
        run(["singularity", "exec", "--bind", str(d), DEMUXAFY_IMG, "Combine_Results.R",
                        "-o", str(d / "old_combined.tsv"), "--demuxlet", str(d / "demuxlet"),
                        "--freemuxlet", str(d / "freemuxlet"), "--souporcell", str(d / "souporcell"),
                        "--vireo", str(d / "vireo"), "--method", "AtLeastHalfSinglet", "--ref", "Vireo",
                        "--pct_agreement", str(pct)], timeout=600)
        if not (d / "old_combined_w_combined_assignments.tsv").exists():
            raise RuntimeError("no combined output")
    except subprocess.CalledProcessError as e:  # Demuxafy itself fails on some corners
        err = (e.stderr or b"").decode(errors="replace")
        if "clusters had at least" in err and "replacement has 0 rows" in err:
            # Known: Combine_Results.R crashes when no cluster of a method reaches
            # pct_agreement. consensus.py treats every such cluster as unmapped.
            return seed, pct, "demuxafy-crash", "no cluster named for a method (known Demuxafy bug)"
        if "allow.cartesian" in err:
            # Known: at pct_agreement <= 0.5 a cluster can map to several donors, and
            # the duplicated join outgrows data.table's cartesian limit.
            return seed, pct, "demuxafy-crash", "duplicate cluster mappings exceed allow.cartesian (known)"
        return seed, pct, "demuxafy-error", err.strip()[-200:]
    rscript(d, "add_vireosnp", old_add_vireosnp,
            'demuxafy_assignments = "%s", vireosnp_assignments = "%s"'
            % (d / "old_combined_w_combined_assignments.tsv", d / "vireosnp/clusters_annotated.tsv"),
            'assignments_with_vireosnp = "%s"' % (d / "old_calls.tsv"))
    wrangle(d, d / "old_calls.tsv", d / "old")
    # new
    sys.path.insert(0, str(REPO / "workflow/scripts"))
    import consensus as c
    t = {k: c.read_table(p) for k, p in dict(
        demuxlet=d / "demuxlet/demuxlet.best", freemuxlet=d / "freemuxlet/freemuxlet.clust1.samples.gz",
        souporcell=d / "souporcell/clusters.tsv", vireo=d / "vireo/donor_ids.tsv",
        vireosnp=d / "vireosnp/clusters_annotated.tsv").items()}
    c.write_tsv(c.build(pct_agreement=pct, **t), d / "new_calls.tsv")
    wrangle(d, d / "new_calls.tsv", d / "new")
    diffs = [f for f in ("sample_barcode_assignments.tsv", "cell_assignment.tsv")
             if (d / "old" / f).read_bytes() != (d / "new" / f).read_bytes()]
    dup = sum(1 for _ in open(d / "old_calls.tsv")) - sum(1 for _ in open(d / "new_calls.tsv"))
    if diffs and apply_bug_applies(d / "old_combined_w_combined_assignments.tsv"):
        # Known: Combine_Results.R takes each barcode's donor from apply(); when every
        # barcode has exactly one distinct donor name, apply() simplifies its result
        # to a vector, the names are lost and every consensus identity becomes NA
        # ("unassigned"). Only reachable in tiny captures; not reproduced.
        return seed, pct, "demuxafy-bug", "apply() simplification: all identities lost"
    return seed, pct, "DIFF " + ",".join(diffs) if diffs else "identical", "dup_rows_collapsed=%d" % dup


def apply_bug_applies(combined):
    import csv
    with open(combined) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    cols = [c for c in rows[0] if c.endswith("_Individual_Assignment") and not c.startswith("AtLeastHalf")]
    return all(len({r[c] for c in cols} - {"doublet", "unassigned"}) == 1 for r in rows)


def safe_run_seed(seed, root, old_add):
    try:
        return run_seed(seed, root, old_add)
    except Exception as e:  # report, don't abort the other seeds
        return seed, "-", "HARNESS-ERROR", repr(e)[-200:]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", default="1-60")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--old-ref", default="origin/main")
    ap.add_argument("--workdir", type=Path)
    a = ap.parse_args()
    lo, hi = map(int, a.seeds.split("-"))
    root = a.workdir or Path(tempfile.mkdtemp())
    root.mkdir(parents=True, exist_ok=True)
    old_add = root / "add_vireosnp_to_combined_results.R"
    old_add.write_bytes(subprocess.run(
        ["git", "-C", str(REPO), "show", a.old_ref + ":workflow/scripts/add_vireosnp_to_combined_results.R"],
        check=True, capture_output=True).stdout)
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        results = sorted(ex.map(lambda s: safe_run_seed(s, root, old_add), range(lo, hi + 1)))
    bad = 0
    for seed, pct, verdict, note in results:
        print("seed %3d  pct=%-5s %-14s %s" % (seed, pct, verdict, note))
        bad += verdict.startswith("DIFF") or verdict in ("HARNESS-ERROR", "demuxafy-error")
    tally = {}
    for _, _, v, _ in results:
        key = v.split()[0]
        tally[key] = tally.get(key, 0) + 1
    print("\n%s   (workdir %s)" % ("  ".join("%s: %d" % kv for kv in sorted(tally.items())), root))
    print("identical = same bytes; demuxafy-crash / demuxafy-bug = known Demuxafy failures the"
          " port does not reproduce; DIFF / HARNESS-ERROR / demuxafy-error = investigate")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
