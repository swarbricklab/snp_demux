# ===================
# title: "count cells"
# description: Counts SNP-based cell clusters (souporcell) with known SNP genotype information (SNP arrays) - singlets as well as all doublets
# project: "workflow-snp-demux"
# author: "Beata Kiedik"
# email: "b.kiedik@garvan.org.au"
# company: "Garvan Institute"
# created: 14/12/2022
# -------------------
# NOTE: this is a snakmake script that uses the "snakemake" object internals
# ===================

import logging as log
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

log.basicConfig(filename = snakemake.log[0], 
                    level = log.INFO, 
                    format = '%(asctime)s:%(levelname)s:%(message)s')                    
#                    format = '%(asctime)s:%(levelname)s:%(name)s:%(message)s')
log.info("LOG FILE: " + snakemake.log[0])

log.info("----------------------------------") 
log.info("CONFIG")

# snakemake driven cmd line parameters
log.info("INPUT")

# souporcell:SNP annotated clusters
plot_clusters_annotated = snakemake.input["clusters_annotated"]
log.info("plot_clusters_annotated: " + plot_clusters_annotated)

log.info("OUTPUT")
plot_cell_counts_all = snakemake.output["cell_counts_all"]
log.info("plot_cell_counts_all: " + plot_cell_counts_all)

plot_cell_counts_singlets = snakemake.output["cell_counts_singlets"]
log.info("plot_cell_counts_singlets: " + plot_cell_counts_singlets)

table_cell_counts = snakemake.output["table_cell_counts_all"]
log.info("table cell counts: " + table_cell_counts)

log.info("")

log.info("WILDCARDS")
capture_id = snakemake.wildcards["capture_id"]
log.info("capture_id: " + capture_id)
log.info("")

log.info("PARAMS")
log.info("")

log.info("----------------------------------") 


# ------------------------------------------------
# Plot number of cells (singlets) before and after annotation per lane 
# ------------------------------------------------

plot_clusters_annotated = pd.read_table(plot_clusters_annotated)

clusters_count = plot_clusters_annotated["MajoritySinglet_Individual_Assignment"].value_counts().rename_axis('sampleID').reset_index(name='cell_counts')
print(clusters_count)
 
#clusters_count = table["MajoritySinglet_Individual_Assignment"].value_counts().rename_axis('sampleID').reset_index(name='cell_counts')

# save tabel with cell counts per sample (with doublets and unassigned)
table_clusters_count = plot_clusters_annotated.groupby(['MajoritySinglet_Individual_Assignment']).size()
table_clusters_count.to_csv(table_cell_counts, sep='\t')

# --------------
# SINGLETS AND DOUBLETS
# --------------

# plot 

def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i, y[i], y[i], ha = 'center')

plt.figure(facecolor='white')

x_pos = np.arange(len(clusters_count["sampleID"]))
 
# Create bars
plt.bar(x_pos, clusters_count["cell_counts"], color=["darkorange"])

# call function to add labels
addlabels(x_pos, clusters_count["cell_counts"])

# Rotation of the bar names
plt.xticks(x_pos, clusters_count["sampleID"], fontsize=5, rotation=45)
#plt.xticks(x_pos, clusters_count["sampleID"], rotation=90)

plt.title('Demuxafy cell count')
plt.xlabel('sampleID / doublet')
plt.ylabel('cell count')

# Show graphic
# plt.show()

# Plot save 
plt.savefig(plot_cell_counts_all)

# --------------
# ONLY SINGLETS
# --------------

# filter only singlets
patternDel = "doublet"
filter = clusters_count["sampleID"].str.contains(patternDel)
clusters_count_singlets = clusters_count[~filter]

def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i, y[i], y[i], ha = 'center')

# plot
plt.figure(facecolor='white')

x_pos = np.arange(len(clusters_count_singlets["sampleID"]))
 
# Create bars
plt.bar(x_pos, clusters_count_singlets["cell_counts"], color=["darkorange"])
 
# call function to add labels
# addlabels(x_pos, clusters_count_singlets["cell_counts"])

# Rotation of the bar names
plt.xticks(x_pos, clusters_count_singlets["sampleID"], fontsize=5, rotation=45)
#plt.xticks(x_pos, clusters_count_singlets["sampleID"], rotation=90)

plt.title('Demuxafy cell count')
plt.xlabel('sampleID')
plt.ylabel('cell count')

# Plot save
plt.savefig(plot_cell_counts_singlets)

# TO DO 
# plot per capture = two lanes facetted into one plot 
# add expected cell count line 
# next plot everything together

log.info("")

# ------------------------------------------------
log.info("PROCESSING COMPLETE")
# ------------------------------------------------

