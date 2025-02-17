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
message("")


message("--------------------------------")
message("COMBINE") # ----
dfCombined <- NULL

for(data_file in data_files)
{
  message("FILE: ", data_file)
  df_data <- read.delim(data_file, stringsAsFactors = F)
  message("Number of rows: ", nrow(df_data))
  
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
