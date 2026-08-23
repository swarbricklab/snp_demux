#!/usr/bin/env Rscript
#' Summarise ambient RNA per barcode using vireo outputs
#' Script is intended to be run via Snakemake.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

# Redirect output to log file
log   <- file(snakemake@log[[1]], open="wt")
sink(log); sink(log, type="message")

# Get file paths from snakemake named inputs/outputs
prop_ambient_file <- snakemake@input[['prop_ambient']]
donor_ids_file    <- snakemake@input[['donor_ids']]
summary_file      <- snakemake@output[['summary']]
stats_file        <- snakemake@output[['stats']]
plot_file         <- snakemake@output[['plot']]

# Read input files
ambient <- read_tsv(prop_ambient_file, show_col_types = FALSE) %>%
    pivot_longer(
        -c("cell", "logLik_ratio"),
        names_to = "donor_id",
        values_to = "prop_ambient"
    ) 
glimpse(ambient)

donors  <- read_tsv(donor_ids_file, show_col_types = FALSE) %>% 
  select(cell, donor_id)
glimpse(donors)

# Filter by left joining on donors to pick out ambient RNA matching vireo assignments
summary <- left_join(
      donors,
      ambient,
      by = c("cell", "donor_id")
    ) %>%
    mutate(
        ambient_other = 1 - prop_ambient
    ) %>%
    select(
      barcode = cell,
      ambient_other,
      logLik_ratio
    )
glimpse(summary)

# Calculate and write summary statistics by binning ambient RNA values
bin_breaks <- c(0.00, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 1.00)
total_barcodes <- summary %>%
    filter(!is.na(ambient_other)) %>%
    nrow()
message("Total barcodes with ambient RNA estimates: ", total_barcodes)

stats <- tibble(
        bin = cut(
            summary$ambient_other,
            breaks = bin_breaks,
            include.lowest = TRUE,
            right = FALSE
        )
    ) %>%
    count(
        bin, name = "count"
    ) %>%
    mutate(
        frac = count / total_barcodes
    ) %>%
    filter(
        !is.na(bin)
    )
glimpse(stats)

# Make bar plot of the stats
#
# Assign the plot rather than leaving it as a bare top-level expression.
# Rscript auto-prints top-level expressions, and printing a ggplot requires a
# graphics device; with none open, R falls back to its default pdf() and writes
# an Rplots.pdf into the working directory -- the dataset workspace root, since
# that is where Snakemake runs jobs. This rule is per-capture, so every capture
# in the run rewrote the same stray file.
p <- ggplot(stats, aes(x = bin, y = frac)) +
  geom_col() +
  labs(
    x = "Ambient RNA fraction",
    y = "Fraction of barcodes",
    title = "Ambient RNA fraction per barcode"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(
    filename = plot_file,
    plot = p,
    width = 8, height = 6, dpi = 300
)

# Write output to CSV
message("Writing summary CSV file")
write_csv(summary, summary_file, na = "")
message("Writing stats CSV file")
write_csv(stats, stats_file, na = "")

