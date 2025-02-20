# ===================
# title: "prepare data for reporting.run after combine_results.R"
# description: extracts cell assignments from combined_results_w_combined_assignments and generates config-input for scRNAseq analysis 
# project: "demuxafy"
# author: "Beata Kiedik"
# email: "b.kiedik@garvan.org.au / b.kiedik@gmail.com"
# company: "Garvan Institute"
# created: 07/03/2023
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
library(tidyr)

message("INPUT")
if(debug) {
    #file_metadata <- "resources/testing/config/samples.mixed_pool12.tsv"
    #file_metadata <- "resources/testing/config/samples.mixed_splitpool.tsv"
    file_metadata <- "resources/testing/config/samples.mixed_pool12.missing_PBMCs.tsv"
} else {
    file_metadata <- snakemake@input[["metadata"]]
}
message("metadata file: ", file_metadata)

if(debug) {
    #file_demuxafy_assignments <- "results/demuxafy/combined_AtLeastHalfSinglet/demuxafy_test_mixed_pool12/Mixed_1/combined_results_w_combined_assignments.inc_vireosnp.tsv"
    #file_demuxafy_assignments <- "results/demuxafy/test_mixed_splitpool/combined_AtLeastHalfSinglet/test_mixed_splitpool/Mixed_1/combined_results_w_combined_assignments.inc_vireosnp.tsv"
    file_demuxafy_assignments <- "results/demuxafy/test_mixed_pool12.missing_PBMCs/combined_AtLeastHalfSinglet/test_mixed_pool12.missing_PBMCs/Mixed_1/combined_results_w_combined_assignments.inc_vireosnp.tsv"    
} else {
    file_demuxafy_assignments <- snakemake@input[["demuxafy_assignments"]]
}
message("file_demuxafy_assignments: ", file_demuxafy_assignments)
message("")

message("OUTPUT")
combined_edited <- snakemake@output[["sample_barcode_assignments"]]
message("combined df: ", combined_edited)
file_all_cells_stats <- snakemake@output[["all_cells_stats"]]
message("all_cells_stats: ", file_all_cells_stats)
file_all_singlets_stats <- snakemake@output[["all_singlets_stats"]]
message("all_singlets_stats: ", file_all_singlets_stats)
file_cell_assignment <- snakemake@output[["cell_assignment"]]
message("file_cell_assignment: ", file_cell_assignment)
message("")

message("WILDCARDS")

if(debug) {
    capture_id <- "Mixed_1"
} else {
    capture_id <- snakemake@wildcards[["capture_id"]]
}
message("capture_id: ", capture_id)
message("")

message("PARAMS")
#MajoritySinglet_Individual_Assignment or AtLeastHalfSinglet_Individual_Assignment
if(debug) {
    assignment_column <- "AtLeastHalfSinglet_Individual_Assignment"
} else {
    assignment_column <- snakemake@params[["assignment_column"]]
}
message("assignment_column: ", assignment_column)

# MajoritySinglet_DropletType or AtLeastHalfSinglet_DropletType
if(debug) {
    droplettype_column <- "AtLeastHalfSinglet_DropletType"
} else {
    droplettype_column <- snakemake@params[["droplettype_column"]]
}
message("droplettype_column: ", droplettype_column)

if(debug) {
    spikein_sample_id <- "PBMC"
} else {
    spikein_sample_id <- snakemake@params[["spikein_sample_id"]]
}
message("spikein_sample_id: ", spikein_sample_id)
message("----------------------------------")

message("Parse meta-data sample-sheet")
df_metadata <- read.delim(file_metadata, header = TRUE, sep = "\t")
#head(df_metadata)

message("read the genotype - cell annotation results")
df_demuxafy_assignments <- read.delim(file_demuxafy_assignments, sep = "\t", header = TRUE)

# get a list of the actual IDs assigned by demuxafy
demuxafy_assignment_ids <- unique(df_demuxafy_assignments$AtLeastHalfSinglet_Individual_Assignment)
demuxafy_assignment_ids <- demuxafy_assignment_ids[!(demuxafy_assignment_ids %in% c("doublet", "unassigned"))]
message("Demuxafy assignment IDs: ", paste0(demuxafy_assignment_ids, ", "))

remap_cell_assignments <- TRUE
cell_assignment_column <- snakemake@config[["cell_assignment_column"]]


message("Identify any samples with no known genotype information")
# samples_without_SNP <- unique(df_metadata$SAMPLE_ID[is.na(df_metadata$GENOTYPE_SAMPLE_ID)])
# warning("Identified ", length(samples_without_SNP), " samples without known genotype information: ", samples_without_SNP)
# samples_without_SNP <- NULL
# samples_without_SNP <- unique(df_metadata[, cell_assignment_column][is.na(df_metadata[, genotype_sample_id_column])])
# warning("Identified ", length(samples_without_SNP), " samples without known genotype information: ", samples_without_SNP)

# should we use exisiting (GENOTYPE_SAMPLE_ID based) assignments or re-map assignments?
# if(cell_assignment_column == genotype_sample_id_column) {
#     message("using exisiting ", genotype_sample_id_column, " based cell assignments...")
#     remap_cell_assignments <- FALSE

# } else {
#     message("re-mapping ", genotype_sample_id_column, " based cell assignments to ", cell_assignment_column, "...")

#     # extract mappings for this capture_id    
#     df_mappings <- df_metadata[df_metadata[, "CAPTURE_ID"] == capture_id, c(cell_assignment_column, genotype_sample_id_column)]
    
#     # check for samples with no known genotype
#     if(!(is.null(samples_without_SNP))) {
#         #no_known_genotype <- df_mappings[is.na(df_mappings[, genotype_sample_id_column]), cell_assignment_column]
#         message("identified samples with no known genotype information, not including them in re-mapping step")
#         df_mappings <- df_mappings[!(is.na(df_mappings[, genotype_sample_id_column])), ]
#     }
    
#     expected_assignment_ids <- df_mappings[, genotype_sample_id_column]
#     message("Expected assignment IDs: ", paste0(expected_assignment_ids, ", "))

#     # extract all mappings
#     #df_mappings <- df_metadata[, c(cell_assignment_column, genotype_sample_id_column)]
    
#     # Identify any (unexpected) cell assignments from samples that are not part of pool
#     not_in_pool <- setdiff(demuxafy_assignment_ids, expected_assignment_ids)
#     message("identified: ", length(not_in_pool), " unexpected cell assignments")
#     if(length(not_in_pool) > 0) {
#         message("these are: ", paste0(not_in_pool, ", "))
#         for(i in not_in_pool) {
#             message("replacing ", i, " with not_in_pool")
#             df_demuxafy_assignments[df_demuxafy_assignments == i] <- "not_in_pool"
#         }
#     }

#     # Now, re-map each of the expected mappings
#     for(i in 1:nrow(df_mappings)) {
#         genotype_sample_id <- df_mappings[i, genotype_sample_id_column]
#         remapped_id <- df_mappings[i, cell_assignment_column]
        
#         message("re-mapping ", genotype_sample_id_column, ": ", genotype_sample_id, " -> ", cell_assignment_column, ": ", remapped_id)
        
#         if(genotype_sample_id == remapped_id) {
#             message(genotype_sample_id_column, ": ", genotype_sample_id, " and ", cell_assignment_column, ": ", remapped_id, " ARE EQUAL - no re-mapping done")

#         } else {
#             message("RE-MAPPING cell annotations from ", genotype_sample_id_column, ": ", genotype_sample_id, " to ", cell_assignment_column, ": ", remapped_id)
#             df_demuxafy_assignments[df_demuxafy_assignments == genotype_sample_id] <- remapped_id
#         }
#     }
# }

# #if(length(samples_without_SNP) > 0) {
# if(!(is.null(samples_without_SNP))) {

#     message("NOTE demuxafy currently doesn't handle sample assignment when no known genotype is provided.")
#     message("So, we default to the assigning souporcell based VireoSNP annotation result to those cells.")

#     # first, record the original demuaxfy annotaion
#     df_demuxafy_assignments$original_demuxafy_Individual_Assignment <- df_demuxafy_assignments[, assignment_column]
#     df_demuxafy_assignments$original_demuxafy_DropletType <- df_demuxafy_assignments[, droplettype_column]

#     mask <- (df_demuxafy_assignments$Souporcell_vireosnp_Individual_Assignment %in% samples_without_SNP)
#     df_demuxafy_assignments[mask, assignment_column] <- df_demuxafy_assignments[mask, "Souporcell_vireosnp_Individual_Assignment"]
#     df_demuxafy_assignments[mask, droplettype_column] <- df_demuxafy_assignments[mask, "Souporcell_vireosnp_DropletType"]    

# } else {
#     message("All samples have known genotype information")
# }

message("write updated demuxafy_assignments to file.")
write.table(df_demuxafy_assignments, file=combined_edited, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

message("------------------------------")
message("stats")
message("------------------------------")

message("Calculate total and fraction of cells assigned to each genotype")
df_assignments <- df_demuxafy_assignments %>% select(paste(assignment_column))
df_assignments_summary <- as.data.frame(table(df_assignments))
  
sum <- sum(df_assignments_summary$Freq)
df_assignments_summary$fraction_all_cells <- paste0(round((df_assignments_summary$Freq)/sum * 100), "%")
colnames(df_assignments_summary) <- c("sample_id", "number_of_cells", "fraction_of_all_cells_captured")

message("Calculate cell singlet totals and as a fraction of cell spike-in content")
# assigned singlets
table_singlets <- df_assignments_summary %>% 
    filter(sample_id != "doublet") %>% 
    filter(sample_id != "singlet") %>% 
    filter(sample_id != "unassigned") %>%
    select(sample_id, number_of_cells)
sum_singlets <- sum(table_singlets$number_of_cells)
table_singlets$fraction_of_all_singlets_captured <- paste0(round((table_singlets$number_of_cells)/sum_singlets * 100), "%")

table_singlets_no_spikein <- table_singlets %>% 
    filter(!(paste(assignment_column) %in% spikein_sample_id)) %>%
    select(sample_id, number_of_cells)
sum_singlets_no_spikein <- sum(table_singlets_no_spikein$number_of_cells)
table_singlets_no_spikein$fraction_of_singlets_captured_without_spikein <- paste0(round((table_singlets_no_spikein$number_of_cells)/sum_singlets_no_spikein * 100), "%")

message("Join the tables")
df_all_singlets <- left_join(table_singlets, table_singlets_no_spikein, by = c("sample_id", "number_of_cells"))
df_all_singlets$capture_id <- paste(capture_id)

df_all_cells <- left_join(df_assignments_summary, df_all_singlets, by = c("sample_id", "number_of_cells"))
df_all_cells$capture_id <- capture_id

message("Extract numbers of non-singlets and spike-ins (if any used)")

# TODO: is this needed? It's not used in the code?
doublets <- df_assignments_summary %>% 
    filter(sample_id == "doublet") %>% 
    pull("number_of_cells")

unassigned <- df_assignments_summary %>% 
    filter(sample_id %in% c("singlet", "unassigned")) %>% 
    pull("number_of_cells")

# TODO: is this needed? It's not used in the code?
spikeins <- df_assignments_summary %>% 
    filter(sample_id == spikein_sample_id) %>% 
    pull("number_of_cells") 

write.table(df_all_cells, file=paste(file_all_cells_stats), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(df_all_singlets, file=paste(file_all_singlets_stats), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

message("------------------------------")
message("Generate scRNAseq workflow config")
message("------------------------------")
result_sc <- df_demuxafy_assignments
scRNAseq_metadata <- result_sc %>% select(Barcode, paste(assignment_column), paste(droplettype_column))

message("change colnames: barcode, assignment (sampleID or doublet), status (singlet, doublet or unassigned)")
colnames(scRNAseq_metadata) <- c("barcode", "assignment", "status")
scRNAseq_metadata$status[grep("singlet", scRNAseq_metadata$assignment)] <- "unassigned"
write.table(scRNAseq_metadata, file=paste(file_cell_assignment), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

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
