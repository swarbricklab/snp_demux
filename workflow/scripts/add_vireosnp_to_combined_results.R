# ===================
# title: "add_vireosnp_to_combined_results"
# description: Adds an extra column to the "demuxafy combined" assignments, based on souporcell::vireosnp annotation results.  
# project: "demuxafy"
# author: "Daniel Roden"
# email: "d.roden@garvan.org.au"
# company: "Garvan Institute"
# created: 22/06/2023
# version: 0.1
# ===================
debug <- FALSE

file_log <- file(snakemake@log[[1]], open="wt")
sink(file_log)
sink(file_log, type="message")

message(paste0("LOG FILE: ", snakemake@log[[1]]))

start_time <- Sys.time()
message("START TIME: ", start_time)
message("\n----------------------------------") 

message("CONFIG")
library(dplyr)

message("INPUT")

if(debug) {
    file_demuxafy_assignments <- "results/demuxafy/combined_AtLeastHalfSinglet/demuxafy_test_mixed_pool12/Mixed_1/combined_results_w_combined_assignments.tsv"
} else {
    file_demuxafy_assignments <- snakemake@input[["demuxafy_assignments"]]
}
message("file_demuxafy_assignments: ", file_demuxafy_assignments)

if(debug) {
    file_vireosnp_assignments <- "results/demuxafy/souporcell_annotation/vireosnp/demuxafy_test_mixed_pool12/Mixed_1/clusters_annotated.tsv"
} else {
    file_vireosnp_assignments <- snakemake@input[["vireosnp_assignments"]]
}
message("file_vireosnp_assignments: ", file_vireosnp_assignments)
message("")

message("OUTPUT")
file_assignments_with_vireosnp <- snakemake@output[["assignments_with_vireosnp"]]
message("file_assignments_with_vireosnp: ", file_assignments_with_vireosnp)
message("")

if(debug) {
    capture_id <- "Mixed_1"
} else {
    capture_id <- snakemake@wildcards[["capture_id"]]
}
message("capture_id: ", capture_id)
message("")

message("PARAMS")
message("")
message("----------------------------------")

message("read the combined results per pool")
df_demuxafy_assignments <- read.delim(file_demuxafy_assignments, sep = "\t", header = TRUE)

# Defensive: Demuxafy's Combine_Results.R can emit exact-duplicate rows (a join
# artifact seen with relaxed correlation/agreement thresholds on sparse data),
# which otherwise propagate into cell_assignment.tsv and break the report's
# per-barcode pivot. These are byte-identical duplicate rows, so dropping them is
# safe. (A barcode with genuinely conflicting rows would NOT be collapsed here.)
n_before <- nrow(df_demuxafy_assignments)
df_demuxafy_assignments <- dplyr::distinct(df_demuxafy_assignments)
if (nrow(df_demuxafy_assignments) < n_before) {
  message("dropped ", n_before - nrow(df_demuxafy_assignments),
          " exact-duplicate demuxafy row(s)")
}

message("read the VireoSNP cell assignments")
df_vireosnp_assignments <- read.delim(file_vireosnp_assignments, sep = "\t", header = TRUE)

# check if demuxafy and vireosnp have same number of cells
# We'd expect these to be the same, unless some additional filtering has taken place
# TODO: In the mixed_pool12 - Mixed_1 results, there is one barcode "TGCGGGTAGACGTCCC-1" in vireosnp but not demuxafy, why, what causes this difference?
barcode_differences <- dplyr::symdiff(df_demuxafy_assignments$Barcode, df_vireosnp_assignments$barcode)
if(length(barcode_differences) > 0) {
    warning("WARNING: The following ", length(barcode_differences), " cell barcodes are present in only one of vireosnp or demuxafy combined: ", paste(barcode_differences))
} else {
    message("All cell barcodes are present in both vireosnp or demuxafy combined")
}

message("Join the tables")
df_with_vireosnp <- left_join(df_demuxafy_assignments, df_vireosnp_assignments[, c("barcode", "status", "assignment")], by = c("Barcode" = "barcode"))
colnames(df_with_vireosnp)[colnames(df_with_vireosnp) %in% c("status", "assignment")] <- c("Souporcell_vireosnp_DropletType", "Souporcell_vireosnp_Individual_Assignment")

message("write updated demuxafy_assignments to file.")
write.table(df_with_vireosnp, file = file_assignments_with_vireosnp, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

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

