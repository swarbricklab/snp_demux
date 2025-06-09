suppressPackageStartupMessages({
    library(tibble)
    library(dplyr)
    library(stringr)
    library(readr)
})


# Get input, output, log and wildcard from snakemake
file_sample_map                <- snakemake@input[["sample_map"]]
file_demuxafy_assignments      <- snakemake@input[["demuxafy_assignments"]]
combined_edited                <- snakemake@output[["sample_barcode_assignments"]]
summary_assignments            <- snakemake@output[["cell_assignment"]]
file_log                       <- snakemake@log[[1]]
capture                        <- snakemake@wildcards[["capture_id"]]
assignment_column              <- snakemake@params[["assignment_column"]]
droplettype_column             <- snakemake@params[["droplettype_column"]]

# Set up logging
log   <- file(snakemake@log[[1]], open="wt")
sink(log); sink(log, type="message")

# Load inputs
df_map     <- read_csv(file_sample_map)
df_assign  <- read_tsv(file_demuxafy_assignments)

# Build donor→sample lookup
map_vec <- df_map %>%
  dplyr::filter(capture_id == capture) %>%
  dplyr::select(donor_id, sample_id) %>%
  distinct() %>%
  deframe()
print(map_vec)

# Recode *_Individual_Assignment columns from donors to samples
# Prefix barcode with capture id
df_assign <- df_assign %>%
  dplyr::mutate(
    across(
      ends_with("_Individual_Assignment"),
      ~ map_vec[.x]
    )
  )

message("Creating assignment summary")
df_summary <- df_assign %>%
  dplyr::select(
    barcode    = Barcode,
    assignment = !!sym(assignment_column),
    status     = !!sym(droplettype_column)
  ) %>%
  dplyr::mutate(
    status = if_else(
        str_detect(assignment, "singlet"), 
        "unassigned", 
        status
    )
  )

message("Writing sample barcode assignments")
write_tsv(df_assign, combined_edited)
message("Writing summary assignments")
write_tsv(df_summary, summary_assignments)
