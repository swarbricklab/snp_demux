# ============================================================================
# negative_test.py  (snakemake per-capture, POST-HOC)
# ----------------------------------------------------------------------------
# Sanity-check a by-elimination donor recovery. A pool's single UNGENOTYPED donor
# X (flagged verdict=ungenotyped by present_absent_screen) is recovered by
# elimination -- it has no genotype, so it can't be positively screened. Here we
# test, using the per-barcode private-alt score matrix M and the final hier
# assignments, whether X's barcodes are a genuinely DISTINCT donor or mislabeled
# known-donor / doublet cells.
#
#   GENOTYPE axis: for each genotyped donor d, pooled private-alt hit-rate over a
#     barcode set S = sum(hit[S,d]) / sum(cov[S,d]). References from THIS pool:
#       own[d]  = rate(cells hier-assigned to d, d)
#       base[d] = rate(cells assigned to OTHER genotyped donors, d)
#       X[d]    = rate(X's barcodes, d)
#       frac[d] = (X[d]-base[d])/(own[d]-base[d])   0=baseline(X!=d) .. 1=own(X==d)
#     frac low for every d -> likely-real; >=HI for one d -> suspect->d';
#     >=HI for >=2 -> doublet/ambient lump.
#   SEX axis (orthogonal; from hier sex_call): X's barcodes should be one sex if X
#     is a single real donor; mixed -> doublet/ambient lump.
#
# The genotype null is NOT zero (donors share private alts at population AF +
# ambient), hence the own/base calibration. Caveat: X is ungenotyped so d's
# "private" sites aren't guaranteed X-negative -> weak for close relatives; the
# sex axis covers that. Emits an n/a row when the pool has no ungenotyped donor.
# ============================================================================
import logging as log
import pandas as pd

log.basicConfig(filename=snakemake.log[0], level=log.INFO,
                format="%(asctime)s:%(levelname)s:%(message)s")

HI = 0.60   # frac >= HI -> X's cells look like that donor (own-level)
LO = 0.35   # frac <= LO -> baseline (X is not that donor)
SEX_CONSISTENT = 0.90

f_pa = snakemake.input["present_absent"]
f_mtx = snakemake.input["score_matrix"]
f_hier = snakemake.input["hier"]
f_out = snakemake.output["tsv"]
capture = snakemake.wildcards["capture_id"]

OUT_COLS = ["capture_id", "ungenotyped", "n_recovered", "genotype_verdict",
            "max_frac_donor", "max_frac", "sex_verdict", "sex_dom_frac"]


def _write(rows):
    pd.DataFrame(rows, columns=OUT_COLS).to_csv(f_out, sep="\t", index=False)
    log.info("wrote %s (%d rows)", f_out, len(rows))


def _na(note):
    log.info("negative_test %s: %s", capture, note)
    # sentinel "none" (NOT "n/a"/"NA" -- those are pandas/R NA tokens on read-back)
    _write([{"capture_id": capture, "ungenotyped": "", "n_recovered": 0,
             "genotype_verdict": "none", "max_frac_donor": "", "max_frac": "",
             "sex_verdict": note, "sex_dom_frac": ""}])


def pooled_rate(M, barcodes, d):
    sub = M.loc[M.index.intersection(barcodes)]
    if sub.empty:
        return None
    cov = sub[f"{d}__cov"].sum()
    return (sub[f"{d}__hit"].sum() / cov) if cov > 0 else None


def main():
    pa = pd.read_csv(f_pa, sep="\t", dtype=str)
    ung = pa.loc[pa["verdict"] == "ungenotyped", "donor"].tolist() if "verdict" in pa else []
    if len(ung) != 1:
        _na("no ungenotyped donor" if not ung else f"multiple ungenotyped {ung}")
        return
    X = ung[0]

    M = pd.read_csv(f_mtx, sep="\t", dtype={"barcode": str})
    if M.empty or "barcode" not in M.columns:
        _na("empty score matrix")
        return
    M = M.set_index("barcode")
    genod = [c[:-5] for c in M.columns if c.endswith("__hit")]

    hier = pd.read_csv(f_hier, sep="\t", dtype={"barcode": str, "donor": str})
    singlet = hier[hier["status"].isin(["singlet", "recovered"])]
    Xbc = singlet.loc[singlet["donor"] == X, "barcode"].tolist()
    n_X = len(Xbc)

    fracs = {}
    for d in genod:
        own = pooled_rate(M, singlet.loc[singlet["donor"] == d, "barcode"].tolist(), d)
        base = pooled_rate(M, singlet.loc[~singlet["donor"].isin([d, X]), "barcode"].tolist(), d)
        xr = pooled_rate(M, Xbc, d)
        if None not in (own, base, xr) and (own - base) > 1e-6:
            fracs[d] = (xr - base) / (own - base)

    high = [d for d, f in fracs.items() if f >= HI]
    max_d = max(fracs, key=fracs.get) if fracs else ""
    max_f = fracs.get(max_d) if max_d else None
    if n_X == 0:
        gv = "no_recovery"
    elif not fracs:
        gv = "untestable"
    elif len(high) >= 2:
        gv = "doublet-lump(" + ",".join(sorted(high)) + ")"
    elif len(high) == 1:
        gv = f"suspect->{high[0]}"
    else:
        gv = "likely-real"

    # sex axis
    sx = hier.loc[hier["barcode"].isin(Xbc), "sex_call"].dropna()
    sx = sx[sx.isin(["male", "female"])]
    sexv, dom = "no_sex", ""
    if len(sx):
        vc = sx.value_counts()
        dom = float(vc.iloc[0] / vc.sum())
        sexv = f"{vc.index[0]}:{dom:.2f}" + ("" if dom >= SEX_CONSISTENT else " MIXED")

    log.info("%s X=%s n=%d gv=%s max=%s(%.3f) sex=%s", capture, X, n_X, gv,
             max_d, (max_f if max_f is not None else -1), sexv)
    _write([{"capture_id": capture, "ungenotyped": X, "n_recovered": n_X,
             "genotype_verdict": gv, "max_frac_donor": max_d,
             "max_frac": round(max_f, 3) if max_f is not None else "",
             "sex_verdict": sexv, "sex_dom_frac": round(dom, 3) if dom != "" else ""}])


if __name__ == "__main__":
    main()
