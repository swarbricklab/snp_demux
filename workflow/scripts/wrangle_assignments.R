suppressPackageStartupMessages({
    library(tibble)
    library(dplyr)
    library(readr)
})

# Demuxafy writes these sentinel values into the *_Individual_Assignment
# columns when it cannot resolve a droplet to a single individual.
# They are not donor IDs, so they must survive the donor -> sample recoding.
demuxafy_sentinels <- c("singlet", "doublet", "unassigned")

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

# Recode a single *_Individual_Assignment column from donors to samples.
# Values that are not donor IDs are passed through unchanged: indexing a named
# vector with an absent name yields NA, so recoding blindly turned every
# sentinel into NA, making "doublet" indistinguishable from a failed lookup.
recode_donors <- function(x) {
  mapped  <- unname(map_vec[x])
  unknown <- !is.na(x) & is.na(mapped) & !(x %in% demuxafy_sentinels)
  if (any(unknown)) {
    warning(
      "Capture ", capture, ": ", sum(unknown),
      " assignments are neither a demuxafy sentinel nor a donor in the sample map: ",
      paste(sort(unique(x[unknown])), collapse = ", ")
    )
  }
  dplyr::coalesce(mapped, x)
}

# Recode *_Individual_Assignment columns from donors to samples
df_assign <- df_assign %>%
  dplyr::mutate(
    across(
      ends_with("_Individual_Assignment"),
      recode_donors
    )
  )

message("Creating assignment summary")
sample_ids <- unname(map_vec)
df_summary <- df_assign %>%
  dplyr::select(
    barcode    = Barcode,
    assignment = !!sym(assignment_column),
    status     = !!sym(droplettype_column)
  ) %>%
  dplyr::mutate(
    # Only real sample IDs belong in the assignment column; the sentinels are
    # reported via 'status'.
    assignment = if_else(assignment %in% sample_ids, assignment, NA_character_),
    # A droplet demuxafy called a singlet but that could not be resolved to a
    # sample in this capture is reported as "unassigned" rather than "singlet".
    status = if_else(
        !is.na(status) & status == "singlet" & is.na(assignment),
        "unassigned",
        status
    )
  )

message("Writing sample barcode assignments")
write_tsv(df_assign, combined_edited)
message("Writing summary assignments")
write_tsv(df_summary, summary_assignments)
