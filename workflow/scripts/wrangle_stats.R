library(dplyr)
library(readr)

# Get input, outputs, params and wildcards from snakemake 
file_assign    <- snakemake@input[["sample_barcode_assignments"]]
file_cells     <- snakemake@output[["all_cells_stats"]]
file_singlets  <- snakemake@output[["all_singlets_stats"]]
assign_col     <- snakemake@params[["assignment_column"]]
dtype_col      <- snakemake@params[["droplettype_column"]]
capture_id     <- snakemake@wildcards[["capture_id"]]

# Load data
df_assign <- read_tsv(file_assign)

# All cells summary. Recovered barcodes (one method's lone singlet call) are
# counted as "recovered", not under their sample, so the per-sample rows (and
# the singlets-only summary) stay consensus singlets.
df_all_cells <- df_assign %>%
  mutate(sample_id = if_else(.data[[dtype_col]] %in% "recovered", "recovered",
                             as.character(.data[[assign_col]]))) %>%
  count(sample_id, name = "number_of_cells") %>%
  mutate(
    fraction_of_all_cells_captured = paste0(
      round(number_of_cells / sum(number_of_cells) * 100), "%"
    ),
    capture_id = capture_id
  )

# Singlets-only summary
df_singlets <- df_all_cells %>%
  filter(!sample_id %in% c("doublet", "unassigned", "singlet", "recovered")) %>%
  mutate(
    fraction_of_all_singlets_captured = paste0(
      round(number_of_cells / sum(number_of_cells) * 100), "%"
    )
  )

# Save results
write_tsv(df_all_cells, file_cells)
write_tsv(df_singlets, file_singlets)
