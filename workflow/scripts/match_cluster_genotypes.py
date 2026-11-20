"""Score every de-novo cluster's genotype against every known donor genotype.

Writes the vireoSNP GPb_diff matrix (mean absolute difference in genotype
probability over shared sites; lower is a closer match) in the same layout as
annotate.py's ``full_GPb_diff_scores_cluster_genotype.tsv``: a ``known genotype``
column, then one ``snp_cluster_<id>`` column per cluster. ``cluster_prefix`` is
stripped from the VCF sample names so ids match the per-barcode calls
(freemuxlet names its VCF samples ``CLUST<id>``; its calls use ``<id>``).

If the two VCFs share no sites the matrix is written with NA scores, so the
consensus leaves every cluster unnamed rather than the workflow failing.
"""
import contextlib
import gzip
import sys

import numpy as np
import pandas as pd
import vireoSNP


def samples(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    return []


def score(known_vcf, cluster_vcf, prefix, log):
    with contextlib.redirect_stdout(log):
        try:
            m = vireoSNP.vcf.match_VCF_samples(known_vcf, cluster_vcf, "GT", "GT")
            donors, clusters, scores = m["full_donors1"], m["full_donors2"], np.asarray(m["full_GPb_diff"])
        except Exception as exc:  # no shared sites, empty VCF, ...
            print("WARNING: genotype matching failed (%r); all scores NA" % exc, file=log)
            donors, clusters = samples(known_vcf), samples(cluster_vcf)
            scores = np.full((len(donors), len(clusters)), np.nan)
    ids = [c[len(prefix):] if prefix and str(c).startswith(prefix) else c for c in clusters]
    df = pd.DataFrame(scores, columns=["snp_cluster_%s" % c for c in ids])
    df.insert(0, "known genotype", [str(d) for d in donors])
    return df


def main(sm):
    with open(sm.log[0], "w") as log:
        df = score(sm.input.known_genotype_vcf, sm.input.vcf, sm.params.cluster_prefix, log)
        df.to_csv(sm.output.scores, sep="\t", index=False, na_rep="NA")
        print(df.to_string(), file=log)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
