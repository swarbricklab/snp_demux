# ===================
# title: "combine_matrices_by_row"
# description: combines all results (row-wise) in input text files
# project: "workflow-single-cell-rna-seq"
# author: "Daniel Roden, Sunny Wu"
# email: "d.roden@garvan.org.au/dlroden@gmail.com, s.wu@garvan.org.au"
# company: "Garvan Institute"
# created: 16/10/2020
# version: 0.1
# ===================
debug <- FALSE

file_log <- file(snakemake@log[[1]], open="wt")
sink(file_log)
sink(file_log, type="message")

message(paste0("LOG FILE: ", snakemake@log[[1]]))
start_time <- Sys.time()
message("START TIME: ", start_time)
message("----------------------------------") 
message("CONFIG")
library(dplyr)
library(data.table)

message("INPUT")
data_files <- snakemake@input[["data_files"]]
message("List of data_files: ", data_files)
message("")

message("OUTPUT")
file_combined <- snakemake@output[["combined"]]
message("file_combined: ", file_combined)
message("")

message("WILDCARDS")
message("")

message("PARAMS")
# Optional. When supplied, one capture_id per element of data_files, in the
# same order -- Snakemake's expand() preserves the order of the capture list,
# so a positional pairing is well defined.
capture_ids <- snakemake@params[["capture_ids"]]
message("capture_ids supplied: ", !is.null(capture_ids))
message("")


message("--------------------------------")
message("COMBINE") # ----
dfCombined <- NULL

if (!is.null(capture_ids) && length(capture_ids) != length(data_files)) {
  stop(
    "capture_ids has ", length(capture_ids), " entries but there are ",
    length(data_files), " input files; they must correspond one-to-one."
  )
}

for(i in seq_along(data_files))
{
  data_file <- data_files[[i]]
  message("FILE: ", data_file)
  df_data <- read.delim(data_file, stringsAsFactors = F)
  message("Number of rows: ", nrow(df_data))

  # Stamp the capture the rows came from. Without this the concatenation is
  # irreversible: sample_barcode_assignments.tsv has no capture column, so the
  # combined table was ~1.8M rows with no way to tell which capture a barcode
  # belonged to, and barcodes are only unique within a capture. Guarded because
  # this rule is also reused for the stats files, which wrangle_stats.R has
  # already stamped.
  if (!is.null(capture_ids) && !("capture_id" %in% colnames(df_data))) {
    df_data$capture_id <- capture_ids[[i]]
  }

  dfCombined <- data.table::rbindlist(list(dfCombined, df_data), fill = TRUE)
}

#dfCombined$sample_id <- stringr::str_split_fixed(dfCombined$cell_ids, pattern = "_", n=Inf)[, 1]
message("number of rows in combined results = ", nrow(dfCombined))
message("writing combined data")
write.table(x = dfCombined, file = file_combined, row.names = FALSE, sep = "\t")

message("--------------------------------")
message("FINALISE") # ----
finish_time <- Sys.time()
print(paste0("start time = ", start_time))
print(paste0("finish time = ", finish_time))

message("--------------------------------")
message("WARNINGS") # ----
warnings()

message("--------------------------------")
message("SESSION INFO") # ----
sessionInfo()
