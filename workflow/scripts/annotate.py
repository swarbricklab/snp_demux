# ===================
# title: "annotation"
# description: Annotates SNP-based cell clusters (souporcell) with known SNP genotype information (SNP arrays)
# project: "workflow-snp-demux"
# author: "Daniel Roden, Nenad Bartonicek"
# email: "d.roden@garvan.org.au"
# company: "Garvan Institute"
# created: 06/06/2022
# -------------------
# NOTE: this is a snakmake script that uses the "snakemake" object internals
# ===================

from distutils.log import error
import numpy as np
import matplotlib.pyplot as plt
import vireoSNP
import pandas as pd
import os
import re
import sys
import logging as log
import itertools

log.basicConfig(filename = snakemake.log[0], 
                    level = log.INFO, 
                    format = '%(asctime)s:%(levelname)s:%(message)s')                    
#                    format = '%(asctime)s:%(levelname)s:%(name)s:%(message)s')
log.info("LOG FILE: " + snakemake.log[0])

# snakemake driven cmd line parameters
log.info("INPUT")

# souporcell derived genotypes of the cell clusters
file_cluster_genotypes = snakemake.input["cluster_genotypes"]
log.info("file_cluster_genotypes: " + file_cluster_genotypes)

# souporcell derived cell -> cluster annotations
file_clusters = snakemake.input["clusters"]
log.info("file_clusters: " + file_clusters)

file_known_genotype_vcf = snakemake.input["known_genotype_vcf"]
log.info("file_known_genotype_vcf: " + file_known_genotype_vcf)
log.info("")

# This is a list of the donors in the pool who have genotype information
# Ultimately this list is derived from the joint VCF file for all donors with genotype information available
genotyped_donors = pd.read_csv(snakemake.input["genotyped_donors"], header=None, dtype=str)[0].tolist()
log.info("Genotyped donors: %s", genotyped_donors)
log.info("")

log.info("OUTPUT")
file_clusters_annotated = snakemake.output["clusters_annotated"]
log.info("file_clusters_annotated: %s", file_clusters_annotated)

file_scores_matched = snakemake.output["scores_matched"]
log.info("file_scores_matched: %s", file_scores_matched)

file_scores_full = snakemake.output["scores_full"]
log.info("file_scores_full: %s", file_scores_full)

file_matching_heatmap_pdf = snakemake.output["matching_heatmap_pdf"]
log.info("file_matching_heatmap_pdf: %s", file_matching_heatmap_pdf)

file_matching_heatmap_png = snakemake.output["matching_heatmap_png"]
log.info("file_matching_heatmap_png: %s", file_matching_heatmap_png)

file_full_heatmap_pdf = snakemake.output["full_heatmap_pdf"]
log.info("file_full_heatmap_pdf: %s", file_full_heatmap_pdf)

file_full_heatmap_png = snakemake.output["full_heatmap_png"]
log.info("file_full_heatmap_png: %s", file_full_heatmap_png)

file_grouped_scores_boxplot = snakemake.output["grouped_scores_boxplot"]
log.info("file_grouped_scores_boxplot: %s", file_grouped_scores_boxplot)

file_final_sample_cluster_mappings = snakemake.output["final_sample_cluster_mappings"]
log.info("file_final_sample_cluster_mappings: %s", file_final_sample_cluster_mappings)

log.info("")

log.info("WILDCARDS")
capture_id = snakemake.wildcards["capture_id"]
log.info("capture_id: %s", capture_id)
log.info("")

log.info("PARAMS")
log.info("")

# internal params
GT_tag0 = 'GT' # common ones: GT, GP, PL
GT_tag1 = 'GT' # common ones: GT, GP, PL

log.info("----------------------------------") 

# Map the samples, via matching of known genotypes to souporcell SNP clusters
vcf_match = vireoSNP.vcf.match_VCF_samples(file_known_genotype_vcf, file_cluster_genotypes, GT_tag0, GT_tag1)
# Donor ids from known genotype VCF that match clusters
matched_donors = vcf_match["matched_donors1"].tolist()
# Souporcell clusters matching these donors
matched_clusters = vcf_match["matched_donors2"].tolist()
# Complete list of souporcell clusters, including clusters not aligned
all_clusters = vcf_match["full_donors2"].tolist()

# ------------------------------------
# MATCHED genotype to cluster mapping
# ------------------------------------
# from here: https://nbviewer.org/github/single-cell-genetics/vireo/blob/master/examples/donor_match.ipynb

log.info("aligned donors:")
log.info(matched_donors)
log.info(matched_clusters)
log.info(all_clusters)

file_config = snakemake.config["deps"]["samples"]
# force the sample IDs to be string to prevent any type errors down-stream (e.g., from numeric sample ids)
df_samples = pd.read_csv(file_config, dtype = {'sample_id': str, 'capture_id': str, 'donor_id': str})
# select the donors associated with the capture pool being processed
donors_in_pool = df_samples.loc[df_samples["capture_id"] == capture_id, "donor_id"].tolist()
log.info("Donors expected in pool: %s", donors_in_pool)

# Default flags for various errors that might occur
unexpected_match = None
contains_unexpected_match = False
unknown_cluster = None
donor_with_unknown_genotype = None
contains_unknown_genotype = False


# Check to see if any matched samples are from another pool (i.e., unexpected) 
unexpected_match = [ donor for donor in matched_donors if donor not in donors_in_pool ]
if len(unexpected_match) > 0:
    contains_unexpected_match = True
    log.error("VireoSNP found matches with samples not expected in the pool: %s", str(unexpected_match))
    raise Exception("ERROR: Match found with sample not in pool")

# Check to see if any samples are without a "known" genotype
donors_with_unknown_genotype = [ donor for donor in donors_in_pool if donor not in genotyped_donors ]

if len(donors_with_unknown_genotype) > 1:
    log.error("No more than one donor in a pool can have an unknown genotype")
    log.error("The following donors do not have a known genotype: %s", donors_with_unknown_genotype)
    raise Exception("ERROR: UNSUPPORTED number of samples without a known genotype in pool")
elif len(donors_with_unknown_genotype) == 1:
    # OK, there is a sample without a "known" genotype
    donor_with_unknown_genotype = donors_with_unknown_genotype[0]
    log.info("Donor %s has no known genotype", donor_with_unknown_genotype)
    contains_unknown_genotype = True
    # NOTE: this isn't needed when using a known genotype VCF file that contains samples accross all captures
    #       However, it is stille required if per-capture known genotype VCFs are used, so leaving in here for legacy support
    # TODO: more testing to ensure the legacy cases still work as expected
    # which SNP cluster is associated with this sample without a genotype?
    unmatched_cluster = [ cluster for cluster in all_clusters if cluster not in matched_clusters ][0]
    unknown_cluster = {unmatched_cluster: donor_with_unknown_genotype}
else:
    log.info("All samples in the pool have a known genotype")

# # Check if we have unexpected matches to samples in other pools and if this coincides with (expected) unknown genotypes in the pool.
if contains_unexpected_match and contains_unknown_genotype:

     # check how many unexpected matches we have
    if len(unexpected_match) > 1:
        log.error("AMBIGUOUS MAPPING: This is an issue as we have clusters assigned to more than one unexepcted sample AND we have a sample with no known genotype information in the pool.")
        log.warning("To allow the workflow to run to completion ALL CLUSTERS mapped to samples not_in_pool have been annotated as the sample with no known genotype information. This is likely incorrect (see below)")
        log.error("This indicates a potential sample drop-out -> try re-running with alternative clustering numbers")
        for i in unexpected_match:
            matched_donors = [donor_with_unknown_genotype if donor == i else donor for donor in matched_donors]
    else:
        log.info("".join(["We have an unexpected match, to sample: ", ",".join(unexpected_match), " and expect a sample: ", ",".join(donor_with_unknown_genotype), ", without a known genotype in the pool"]))
        log.info("assume that this is the missing sample and update the vireoSNP.vcf.match_VCF_samples annotation accordingly")
        matched_donors = [donor_with_unknown_genotype if donor == unexpected_match else donor for donor in matched_donors]

elif contains_unexpected_match and not contains_unknown_genotype:
    log.warning("Found an unexpected match, but there aren't any missing genotypes => this indicates a possible issue (e.g., sample-swap/drop-out)")
elif contains_unknown_genotype and not contains_unexpected_match:
    # Do nothing
    log.info("We have an unknown genotype but no match to another sample. This is OK and is expected if we only use known genotypes from the sample pool (i.e., not integrated across all captures/pools)")
else:
    # This is the majority case, where all samples have known genotypes and all matches are from the known genotypes -> Do nothing
    log.info("All samples in the pool have a known genotype")

# Pull together matched SNP clusters and matched donors
# NOTE: the changed order of donors 1 and 2
dict_matching_donors = dict(zip(matched_clusters, matched_donors))

# update unknown genotype clusters if needed
# NOTE (as above): this isn't needed when using a known genotype VCF file that contains samples accross all captures
#       However, it is still required if per-capture known genotype VCFs are used, so leaving in here for legacy support
if unknown_cluster is not None:
    dict_matching_donors.update(unknown_cluster)

# ------------------------------------------
# Extract the GPb_diff scores to file (i.e., the distance between cluster and known genotype)
# ------------------------------------------
# Matched samples score matrix
df_matched_scores = pd.DataFrame(vcf_match["matched_GPb_diff"])
df_matched_scores.columns = ["snp_cluster_" + x for x in matched_clusters]
df_matched_scores.insert(0, "known genotype", matched_donors)
df_matched_scores.to_csv(file_scores_matched, sep='\t', index=False, header=True)

# Full score matrix
df_full_scores = pd.DataFrame(vcf_match["full_GPb_diff"])
df_full_scores.columns = ["snp_cluster_" + x for x in vcf_match["full_donors2"]]
df_full_scores.insert(0, "known genotype", vcf_match["full_donors1"])
df_full_scores.to_csv(file_scores_full, sep='\t', index=False, header=True)

# keep an unaltered copy of this for later
df_full_scores_orig = df_full_scores.copy()

# -----------------------
# Plot heatmaps
# -----------------------
# Matched samples score matrix
fig = plt.figure()
fig.set_size_inches(10, 8)
vireoSNP.plot.heat_matrix(vcf_match["matched_GPb_diff"], matched_donors, matched_clusters)
plt.title("Geno Prob Delta: MATCHED known genotype vs SNP clusters")
plt.xlabel("Souporcell SNP cluster IDs")
plt.ylabel("Known Genotypes")
plt.tight_layout()
plt.savefig(file_matching_heatmap_pdf)
plt.savefig(file_matching_heatmap_png)

# Full score matrix
fig = plt.figure()
fig.set_size_inches(10, 8)
vireoSNP.plot.heat_matrix(vcf_match["full_GPb_diff"], vcf_match["full_donors1"], vcf_match["full_donors2"])
plt.title("Geno Prob Delta: FULL known genotype vs SNP clusters")
plt.xlabel("Souporcell SNP cluster IDs")
plt.ylabel("Known Genotypes")
plt.tight_layout()
plt.savefig(file_full_heatmap_pdf)
plt.savefig(file_full_heatmap_png)

# ------------------------------------------
# Have any samples dropped out from the capture pool? 
# ------------------------------------------
# NOTE: This should result in:
#  A) One (or more) SNP clusters mapping to a known genotype with poor (i.e., high) scores, 
#     BUT with no strong match to any other known genotype. 
#     This usually indicates that the sample has "mostly" dropped out and the number of cells in the SNP cluster is very low and therefore causing a low probability mapping (high) score
#  B) one (or more) SNP clusters mapping to a known genotype with poor (i.e., high) scores
#     AND with a strong match to another known genotype. 
#     This usually indicates that a sample has dropped out and the SNP cluster is actually from another known genotype and therefore causing a high probability mapping (low) score to an already mapped known genotype (i.e., a "known" SNP cluster has been split in two)
# ------------------------------------------
# NOTE: Let's use a GPb score difference of 0.1 as a (bit of an arbritrary) threshold for now
# TODO: add as an input parameter via config?
score_threshold = 0.1

# first let's get the avg score for "matched" and "un-matched" SNP clusters to the expected known genotypes
matched_scores = []
matched_scores_pairs = []
unmatched_scores = []

for cluster, known in dict_matching_donors.items():
    matched_scores.extend(list(df_full_scores.loc[df_full_scores["known genotype"] == known, "snp_cluster_" + str(cluster), ]))
    matched_scores_pairs.append([known, cluster])
    df_full_scores.loc[df_full_scores["known genotype"] == known, "snp_cluster_" + str(cluster)] = np.nan

df_full_scores = df_full_scores.set_index("known genotype")
unmatched_scores = df_full_scores.to_numpy(np.float64).flatten()
unmatched_scores = unmatched_scores[np.isnan(unmatched_scores) == False]

median_matched = np.median(matched_scores)
median_unmatched = np.median(unmatched_scores)

log.info("MEDIAN MATCH SCORE     = " + str(median_matched))
log.info("MEDIAN UNMATCHED SCORE = " + str(median_unmatched))
log.info("DIFFERENCE IN SCORES   = " + str(median_unmatched - median_matched))

# OK, let's plot these scores
df_grouped_scores = pd.DataFrame({"group": "matched", "score": matched_scores})
df_grouped_scores = pd.concat([df_grouped_scores, pd.DataFrame({"group": "unmatched", "score": unmatched_scores})])
fig = plt.figure()
#fig.set_size_inches(10, 8)
bp = df_grouped_scores.groupby(['group']).boxplot()
plt.savefig(file_grouped_scores_boxplot)

# --------------
# OUTCOME A
# --------------
# Are there any outlier "matching" scores within a specific threshold of the "unmatched" average score?
matched_scores = np.asarray(matched_scores)
matched_outliers = None
matched_outliers_mask = matched_scores > median_unmatched - score_threshold
matched_outliers = matched_scores[matched_outliers_mask]

if len(matched_outliers) == 1:
    # Identify which cluster (and matching genotype) this outlier is associated with
    matched_outliers_pairs = list(*itertools.compress(matched_scores_pairs, matched_outliers_mask))
    snp_cluster = str(matched_outliers_pairs[1])
    known_genotype = str(matched_outliers_pairs[0])

    log.info("MATCHED OUTLIER FOUND: This indicates a poor match between known genotype and SNP cluster")
    log.info("known genotype=" + known_genotype + "; SNP-cluster=" + snp_cluster + "; score=" + str(*matched_outliers))

    # --------------
    # OUTCOME B
    # --------------
    # Does this SNP cluster have a strong matche with any other known genotype?  
    # (i.e., are there any outlier "unmatched" scores within a specific threshold of the "matched" average score?)
    unmatched_outliers_mask = df_full_scores.loc[:, "snp_cluster_" + snp_cluster] < median_unmatched - score_threshold

    if sum(unmatched_outliers_mask) == 1:
        # Identify which matching genotype this un-matched SNP cluster outlier is associated with
        known_genotype = str(unmatched_outliers_mask[unmatched_outliers_mask].index[0])

        # -------------------------------------
        # DEBUG
        # log.info("DEBUG: df_full_scores: ")
        # log.info(df_full_scores)
        # log.info("DEBUG: : type(known_genotype)")
        # log.info(type(known_genotype))
        # log.info("DEBUG: type(df_full_scores.index)")
        # log.info(type(df_full_scores.index))
        # log.info("DEBUG: type(df_full_scores.index[0])")
        # log.info(type(df_full_scores.index[0]))
        # log.info("")
        # -------------------------------------

        matched_score = df_full_scores.loc[known_genotype, "snp_cluster_" + snp_cluster]
        log.info("UN-MATCHED OUTLIER FOUND: This indicates a good match between a known genotype and 2 SNP clusters")
        log.info("known genotype=" + known_genotype + "; SNP-cluster=" + snp_cluster + "; score=" + str(matched_score))
        log.info("UPDATING SNP-cluster to known genotype mapping")

        # OK, let's update the SNP-cluster to known genotype annotations
        log.warning("SAMPLE " + str(dict_matching_donors[snp_cluster]) + " looks like it has dropped out of the capture pool")
        log.warning("REPLACING it with sample " + str(known_genotype))
        dict_matching_donors[snp_cluster] = known_genotype

    elif sum(unmatched_outliers_mask) == 0:
        # OK, there are no better matches to this SNP-cluster,so just keep as-is, but give a warning
        log.warning("There are no better matches to this SNP-cluster, so mapping will remain unchanged")
        log.warning("BUT, it looks like SAMPLE " + str(dict_matching_donors[snp_cluster]) + " has dropped out of the capture pool")

    else:
        log.error("MORE THAN ONE ALTERNATIVE GENOTYPE MATCH FOUND. This means that the SNP cluster mapping to known genotypes is ambiguous. Interpret results of Demux workflow WITH CAUTION.")
        raise Exception("ERROR: MORE THAN ONE ALTERNATIVE GENOTYPE MATCH FOUND. This means that the SNP cluster mapping to known genotypes is ambiguous. Interpret results of Demux workflow WITH CAUTION.")

elif len(matched_outliers) > 1:
    # Identify which cluster (and matching genotype) this outlier is associated with
    matched_outliers_pairs = list(itertools.compress(matched_scores_pairs, matched_outliers_mask))    
    log.error("MORE THAN ONE MATCHED OUTLIER FOUND. This means that more than one SNP cluster is not mapping to a known genotype. Interpret results of Demux workflow WITH CAUTION.")
    log.error(str(matched_outliers_pairs))
    log.error("This indicates a potential sample drop-out -> try re-running with alternative clustering numbers")
    #raise Exception("ERROR: MORE THAN ONE MATCHED OUTLIER FOUND. This means that more than one SNP cluster is not mapping to a known genotype. Demux workflow can't continue.")

elif len(matched_outliers) == 0:
    log.info("NO MATCHED OUTLIERS FOUND -> ALL SNP clusters and known genotypes matched")
    # NOTE: no change to the matching donors dictionary

else:
    log.error("UNEXPECTED number of matched outliers found")
    raise Exception("ERROR: UNEXPECTED number of matched outliers found : " + str(len(matched_outliers)))

# ------------------------------------------------
# Annotate the souporcell SNP-based cluster file
# ------------------------------------------------
df_clusters = pd.read_table(file_clusters)
df_annotated = df_clusters.replace({"assignment": dict_matching_donors})
df_annotated.to_csv(file_clusters_annotated, sep='\t', index=False, header=True)

log.info("")
log.info("FINAL SNP CLUSTER to KNOWN GENOTYPE SAMPLE MAPPINGS")
log.info("---------------------------------------------------")
df_matched_final = pd.DataFrame()
for snp_cluster, sample_id in dict_matching_donors.items():
    gpb_score = df_full_scores_orig.loc[df_full_scores_orig["known genotype"] == sample_id, "snp_cluster_" + str(snp_cluster)]
    
    # check if no match is available (e.g., due to one of the samples not having any known genotype information available)
    if gpb_score.empty:
        gpb_score = np.NaN
    else:
        gpb_score = float(gpb_score)

    log.info("SNP CLUSTER: " + str(snp_cluster) + " <==> " + "SAMPLE ID: " + str(sample_id) + " (GPb_diff_score = " + str(gpb_score) + ")")
    df_row = pd.DataFrame({'SAMPLE_ID': [str(sample_id)], 'SNP_CLUSTER' : [str(snp_cluster)], 'GPb_diff_score' : [gpb_score]})
    df_matched_final = pd.concat([df_matched_final, df_row])

log.info("")

# write to file
df_matched_final.to_csv(file_final_sample_cluster_mappings, sep='\t', index=False, header=True)

# ------------------------------------------------
log.info("PROCESSING COMPLETE")
# ------------------------------------------------
