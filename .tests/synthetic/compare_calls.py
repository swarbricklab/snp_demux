#!/usr/bin/env python3
"""Compare per-method and consensus calls between two snp_demux result dirs
(e.g. golden vs a PR run that is expected to DIFFER). Prints, per capture,
agreement of (DropletType, Individual_Assignment) per method and of the
consensus cell_assignment.tsv (status, assignment), plus a confusion summary.

Usage: compare_calls.py GOLDEN_DIR FRESH_DIR [capture]"""
import sys
import pandas as pd

g, f = sys.argv[1], sys.argv[2]
cap = sys.argv[3] if len(sys.argv) > 3 else "onek1k"
rd = lambda p: pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
sg = rd(f"{g}/assignment/{cap}/sample_barcode_assignments.tsv").drop_duplicates("Barcode").set_index("Barcode")
sf = rd(f"{f}/assignment/{cap}/sample_barcode_assignments.tsv").drop_duplicates("Barcode").set_index("Barcode")
common = sg.index.intersection(sf.index)
print(f"barcodes: golden={len(sg)} fresh={len(sf)} common={len(common)}")
for m in ["Demuxlet", "Freemuxlet", "Souporcell", "Vireo", "AtLeastHalfSinglet", "Souporcell_vireosnp"]:
    a = sg.loc[common, [f"{m}_DropletType", f"{m}_Individual_Assignment"]].agg("|".join, axis=1)
    b = sf.loc[common, [f"{m}_DropletType", f"{m}_Individual_Assignment"]].agg("|".join, axis=1)
    print(f"  {m:22s} agreement {(a == b).mean():7.2%}  ({(a != b).sum()} changed)")
cg = rd(f"{g}/assignment/{cap}/cell_assignment.tsv").set_index("barcode")
cf = rd(f"{f}/assignment/{cap}/cell_assignment.tsv").set_index("barcode")
k = cg.index.intersection(cf.index)
a = cg.loc[k, "status"] + "|" + cg.loc[k, "assignment"]
b = cf.loc[k, "status"] + "|" + cf.loc[k, "assignment"]
print(f"consensus cell_assignment (status|assignment) agreement {(a == b).mean():.2%} ({(a != b).sum()} changed)")
ch = pd.DataFrame({"golden": cg.loc[k, "status"], "fresh": cf.loc[k, "status"]})[a != b]
print("  changed, by status golden->fresh:", ch.value_counts().to_dict())
sing = (cg.loc[k, "status"] == "singlet") & (cf.loc[k, "status"] == "singlet")
swap = (cg.loc[k, "assignment"] != cf.loc[k, "assignment"]) & sing
print(f"  singlet->different sample (donor swaps): {int(swap.sum())}")
print("  status counts golden:", cg["status"].value_counts().to_dict())
print("  status counts fresh: ", cf["status"].value_counts().to_dict())
