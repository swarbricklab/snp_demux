#!/usr/bin/env python3
"""Check a missing-genotype run (scenario missing1) against the golden.

The golden assigns every cell with every donor genotyped. With one donor's
genotype removed, the workflow must still assign that donor by elimination, and
the rest of the capture should barely move. Checked on cell_assignment.tsv:

  * recall of each dropped donor: of the golden singlets assigned to it, the
    fraction the scenario also assigns to it          >= MIN_RECALL (0.80)
  * overall agreement on (status, assignment) over the barcodes both runs
    report                                            >= MIN_AGREEMENT (0.90)

Thresholds are initial guesses for the shallow OneK1K fixture; override with the
environment variables of the same name once calibrated.

Usage: compare_missing.py GOLDEN_CELL_ASSIGNMENT FRESH_CELL_ASSIGNMENT DONOR[,DONOR]
"""
import csv
import os
import sys


def read(path):
    with open(path, newline="") as fh:
        return {r["barcode"]: (r["status"], r["assignment"]) for r in csv.DictReader(fh, delimiter="\t")}


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        return 2
    golden, fresh = read(sys.argv[1]), read(sys.argv[2])
    dropped = sys.argv[3].split(",")
    min_recall = float(os.environ.get("MIN_RECALL", 0.80))
    min_agree = float(os.environ.get("MIN_AGREEMENT", 0.90))

    shared = golden.keys() & fresh.keys()
    agree = sum(golden[b] == fresh[b] for b in shared) / len(shared) if shared else 0.0
    ok = agree >= min_agree
    print(f"barcodes: golden {len(golden)}  fresh {len(fresh)}  shared {len(shared)}")
    print(f"overall agreement: {agree:.4f} (min {min_agree})  {'OK' if agree >= min_agree else 'LOW'}")

    for d in dropped:
        g = [b for b in shared if golden[b] == ("singlet", d)]
        hit = sum(fresh[b] == ("singlet", d) for b in g)
        recall = hit / len(g) if g else 0.0
        n_fresh = sum(v == ("singlet", d) for v in fresh.values())
        good = bool(g) and recall >= min_recall
        ok &= good
        print(
            f"dropped donor {d}: golden singlets {len(g)}  recovered {hit}  "
            f"recall {recall:.4f} (min {min_recall})  scenario singlets {n_fresh}  "
            f"{'OK' if good else 'LOW'}"
        )

    print("\nPASS" if ok else "\nFAIL: missing-genotype run diverges from the golden")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
