#!/usr/bin/env python3
"""Run the doublet-stack scripts (#59 posteriors, #60/#61 summary) on simulated
captures with a known doublet rate, and check the estimates against the truth.

For each replicate seed: generate captures with make_doublet_captures.py, run
doublet_posteriors.py then doublet_summary.py (per capture + cohort roll-up)
through their real `main(snakemake)` entry points with a stand-in snakemake
object, and collect the summary rows.

Reports per scenario, over replicates:
  * f_x_hard mean vs the expected cross-donor fraction  d * (1 - sum p^2) * sens
  * d_total mean vs the true d (unbiased if sensitivity = 1)
  * coverage of d_total +/- 1.96 * d_total_sd (should be ~95% when sens = 1)
  * fallbacks (single-sample capture must give NA)
Plus a single legacy-format run (status NA for non-singlets) to show what a
pre-#22 cell_assignment.tsv does to the stack.

Usage: doublet_stack_check.py WORKTREE [--reps 100]
  WORKTREE is a checkout of the PR branch (needs workflow/scripts/doublet_*.py).
"""
import argparse
import csv
import importlib.util
import io
import math
import statistics
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_main(mod, **kw):
    err = sys.stderr
    try:
        mod.main(SimpleNamespace(**kw))
    finally:
        if sys.stderr is not err:
            sys.stderr.close()
        sys.stderr = err


def run_stack(post, summ, data, rate, ca_name="cell_assignment.tsv"):
    caps = [r["capture_id"] for r in csv.DictReader(open(data / "truth.tsv"), delimiter="\t")]
    outs = []
    for cap in caps:
        c = data / cap
        run_main(post, input=SimpleNamespace(samples=str(data / "samples.csv"), cell_assignment=str(c / ca_name),
                                             vireo=str(c / "donor_ids.tsv")),
                 output=SimpleNamespace(csv=str(c / "doublet_posteriors.csv")), log=[str(c / "post.log")],
                 wildcards=SimpleNamespace(capture_id=cap), rule="doublet_posteriors")
        kw = dict(input=SimpleNamespace(samples=str(data / "samples.csv"), posteriors=str(c / "doublet_posteriors.csv")),
                  output=SimpleNamespace(csv=str(c / "doublet_summary.csv")), log=[str(c / "summ.log")],
                  wildcards=SimpleNamespace(capture_id=cap), rule="doublet_summary")
        if rate is not None:
            kw["params"] = SimpleNamespace(rate_per_1k=rate)
        run_main(summ, **kw)
        outs.append(str(c / "doublet_summary.csv"))
    run_main(summ, input=SimpleNamespace(summaries=outs), output=SimpleNamespace(csv=str(data / "cohort.csv")),
             log=[str(data / "combine.log")], wildcards=SimpleNamespace(), rule="combine_doublet_summary")
    return {r["capture_id"]: r for r in csv.DictReader(open(data / "cohort.csv"))}


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("worktree", type=Path)
    ap.add_argument("--reps", type=int, default=100)
    a = ap.parse_args()
    scripts = a.worktree / "workflow/scripts"
    post = load(scripts / "doublet_posteriors.py", "doublet_posteriors")
    summ = load(scripts / "doublet_summary.py", "doublet_summary")
    has_rates = "expected_rates" in dir(summ)
    rate = 0.008 if has_rates else None

    results, truth = {}, None
    with tempfile.TemporaryDirectory() as tmp:
        for seed in range(1, a.reps + 1):
            d = Path(tmp) / f"s{seed}"
            subprocess.run([sys.executable, str(HERE / "make_doublet_captures.py"), str(d), "--seed", str(seed)], check=True)
            if truth is None:
                truth = {r["capture_id"]: r for r in csv.DictReader(open(d / "truth.tsv"), delimiter="\t")}
            for cap, row in run_stack(post, summ, d, rate).items():
                results.setdefault(cap, []).append(row)
            if seed == 1:
                first = run_stack(post, summ, d, rate)
                stat_vals = {cap: sorted({r["consensus_status"] for r in csv.DictReader(open(d / cap / "doublet_posteriors.csv"))})
                             for cap in truth}
                legacy = run_stack(post, summ, d, rate, ca_name="cell_assignment.legacy.tsv")
                legacy_vals = sorted({r["consensus_status"] for r in csv.DictReader(open(d / "equal4_6pct" / "doublet_posteriors.csv"))})
            shutil.rmtree(d)  # keep the temp footprint to one replicate

    print(f"replicates: {a.reps}   expected-rate columns present: {has_rates}\n")
    hdr = f"{'scenario':12s} {'n_samp':>6s} {'sum_p_sq(true/obs)':>19s} {'f_x exp':>8s} {'f_x obs':>8s}"
    if has_rates:
        hdr += f" {'d_true':>6s} {'d_total':>8s} {'bias%':>6s} {'cover95':>7s} {'fallback':>15s}"
    print(hdr)
    ok = True
    for cap, rows in results.items():
        t = truth[cap]
        fx = statistics.mean(f(r["f_x_hard"]) for r in rows)
        sp = statistics.mean(f(r["sum_p_sq"]) for r in rows)
        line = f"{cap:12s} {rows[0]['n_samples']:>6s} {float(t['sum_p_sq_true']):>9.4f}/{sp:<9.4f} {float(t['f_x_expected']):>8.4f} {fx:>8.4f}"
        if has_rates:
            dt = [f(r["d_total"]) for r in rows]
            fb = {r["d_fallback"] for r in rows}
            d_true = float(t["d_true"])
            if all(math.isnan(x) for x in dt):
                line += f" {d_true:>6.3f} {'NA':>8s} {'':>6s} {'':>7s} {','.join(sorted(fb)):>15s}"
            else:
                cover = statistics.mean(abs(f(r["d_total"]) - d_true) <= 1.96 * f(r["d_total_sd"]) for r in rows)
                bias = 100 * (statistics.mean(dt) - d_true) / d_true
                line += f" {d_true:>6.3f} {statistics.mean(dt):>8.4f} {bias:>+6.1f} {cover:>7.0%} {','.join(sorted(fb)):>15s}"
                if float(t["sensitivity"]) == 1 and float(t["unassigned"]) == 0 and (abs(bias) > 3 or cover < 0.88):
                    ok = False
                    line += "  <-- CHECK"
        print(line)

    print("\nconsensus_status values in doublet_posteriors.csv (seed 1):",
          {k: v for k, v in stat_vals.items() if k in ("equal4_6pct",)})
    print("\nLegacy cell_assignment.tsv (status NA for non-singlets), seed 1, equal4_6pct; consensus_status values:", legacy_vals)
    for col in ("n_doublet_hard", "n_unassigned", "f_x_hard", "d_total", "d_fallback") if has_rates else ("n_doublet_hard", "n_unassigned", "f_x_hard"):
        print(f"   {col:15s} current={first['equal4_6pct'][col]:>10s}   legacy={legacy['equal4_6pct'][col]:>10s}")
    print("\nOVERALL:", "PASS" if ok else "CHECK FLAGGED ROWS")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
