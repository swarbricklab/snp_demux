# ============================================================================
# aggregate_negative_test.py  (snakemake, cohort-level PASS 2)
# ----------------------------------------------------------------------------
# Concatenate the per-capture negative-test rows into one cohort table, then
# FINALISE the relatedness verdict across co-pooled captures.
#
# Per-capture negative_test.py emits a provisional related/suspect call for each
# single-donor "suspect" (frac>=HI for one donor D), split by the LOCAL het_share
# at D's hom-alt sites. But relatedness is a property of the donor PAIR (X, D), not
# of any one capture, and souporcell's de-novo cluster consensus collapses het->hom
# when X has few cells -- so low coverage only ever DEFLATES het_share, never
# inflates it. The least-biased estimate for a pair is therefore the MAX het_share
# across every capture that recovered that pair. Here we pool by (ungenotyped, D),
# take the best-powered het_share, and rewrite each verdict accordingly. Adds a
# `pair_het_share` column recording the value the call was made on.
# ============================================================================
import logging as log
import re
import pandas as pd

log.basicConfig(filename=snakemake.log[0], level=log.INFO,
                format="%(asctime)s:%(levelname)s:%(message)s")

REL_HET = 0.15   # keep in sync with negative_test.py

inputs = list(snakemake.input)
f_summary = snakemake.output["summary"]
f_ckpt = snakemake.output["checkpoint_file"]

_VERDICT_RE = re.compile(r"^(related|suspect)->(.+)$")


def main():
    df = pd.concat([pd.read_csv(f, sep="\t", dtype=str) for f in inputs],
                   ignore_index=True)
    df["pair_het_share"] = ""

    # rows carrying a provisional related/suspect verdict + a parseable candidate D
    def parse_D(v):
        m = _VERDICT_RE.match(str(v))
        return m.group(2) if m else None

    df["_D"] = df["genotype_verdict"].map(parse_D)
    cand = df[df["_D"].notna()].copy()
    cand["_hs"] = pd.to_numeric(cand["het_share"], errors="coerce")

    # PASS 2: best-powered (max) het_share per (ungenotyped X, candidate D) pair
    pair_best = (cand.dropna(subset=["_hs"])
                     .groupby(["ungenotyped", "_D"])["_hs"].max().to_dict())

    for i, row in cand.iterrows():
        best = pair_best.get((row["ungenotyped"], row["_D"]))
        D = row["_D"]
        if best is None:
            # no capture of this pair had a usable het_share -> stay conservative
            verdict = f"suspect->{D}"
            df.at[i, "pair_het_share"] = ""
        else:
            verdict = (f"related->{D}" if best >= REL_HET else f"suspect->{D}")
            df.at[i, "pair_het_share"] = round(best, 3)
        df.at[i, "genotype_verdict"] = verdict
        log.info("finalise %s X=%s D=%s local_hs=%s pair_hs=%s -> %s",
                 row["capture_id"], row["ungenotyped"], D,
                 row.get("het_share"), df.at[i, "pair_het_share"], verdict)

    df = df.drop(columns=["_D"])
    df.to_csv(f_summary, sep="\t", index=False)
    log.info("wrote %s (%d rows)", f_summary, len(df))
    with open(f_ckpt, "w") as fh:
        fh.write("\n".join(inputs) + "\n")


if __name__ == "__main__":
    main()
