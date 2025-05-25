library(dplyr)
library(readr)

# Get input, outputs, params and wildcards from snakemake 
file_assign    <- snakemake@input[["sample_barcode_assignments"]]
file_cells     <- snakemake@output[["all_cells_stats"]]
file_singlets  <- snakemake@output[["all_singlets_stats"]]
assign_col     <- snakemake@params[["assignment_column"]]
capture_id     <- snakemake@wildcards[["capture_id"]]

# Load data
df_assign <- read_tsv(file_assign)

# All cells summary
df_all_cells <- df_assign %>%
  count(!!sym(assign_col), name = "number_of_cells") %>%
  rename(sample_id = !!sym(assign_col)) %>%
  mutate(
    fraction_of_all_cells_captured = paste0(
      round(number_of_cells / sum(number_of_cells) * 100), "%"
    ),
    capture_id = capture_id
  )

# Singlets-only summary
df_singlets <- df_all_cells %>%
  filter(!sample_id %in% c("doublet", "unassigned", "singlet")) %>%
  mutate(
    fraction_of_all_singlets_captured = paste0(
      round(number_of_cells / sum(number_of_cells) * 100), "%"
    )
  )

# Save results
write_tsv(df_all_cells, file_cells)
write_tsv(df_singlets, file_singlets)