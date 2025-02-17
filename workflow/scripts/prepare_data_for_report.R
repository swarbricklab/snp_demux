# ===================
# title: "prepare data for reporting.run after combine_results.R"
# description: extracts cell assignments from combined_results_w_combined_assignments, changes names from CEL based to SampleIDs and generates config-input for scRNAseq analysis 
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

message("INPUT")
metadata <- snakemake@input[["metadata"]]

message("metadata file: ", metadata)

message("OUTPUT")
df_all <- snakemake@output[["df_all"]]
message("combined df: ", df_all)
df_singlets <- snakemake@output[["df_singlets"]]
message("combined df of singlets: ", df_singlets)

message("PARAMS")
threshold <- snakemake@params[["threshold"]] 
#MajoritySinglet_Individual_Assignment or AtLeastHalfSinglet_Individual_Assignment
droplet <- snakemake@params[["droplet"]] 
# MajoritySinglet_DropletType or AtLeastHalfSinglet_DropletType
path <- snakemake@params[["path"]]
message("path: ", path)

message("----------------------------------")
message("EXTRACT SINGLET CELLs in SAMPLEs") # ----

library(dplyr)
library(tidyr) 

metadata <- read.table(metadata, header = TRUE)

# edit metadata to extract the CEL based sample names in order to match them to samples of origin 
# TODO: I think use label based indexing, rather than index/number-based, as it's easy for the column orders to be changed and introduce subtle bugs.
# TODO: Regardless, all of this should be abstracted out of this script anyway
CEL <- metadata[,7]
# some CEL files dont have the same name pattern = between UKB and Results are is no year, so I added it to be able to extract the sample name in the same way for all CEL files
CEL <- gsub("UKB_RESULTS", "UKB_YYYY_RESULTS", CEL)
CEL <- sub('.*^([^_]+_[^_]+_[^_]+_[^_]+_[^_]+_[^_]+_[^_]+_)', '\\2', CEL)
CEL <- sub("\\..*", "", CEL)
CEL <- gsub("5%_PBMC_spike_ins", "PBMC-Donor-6", CEL)
CEL <- stringr::str_replace_all(CEL, "_", "-")
x <- is.na(CEL) 
CEL[x] <- metadata[x,"SAMPLE_ID"]

# add CEL name column to the metadata
metadata$CEL_sample <- CEL

# change PBMC name to a name without %
metadata$SAMPLE_ID <- gsub("5%_PBMC_spike_ins", "PBMC_spike", metadata$SAMPLE_ID)  

capture_ids <- unique(metadata$CAPTURE_ID)

# read the combined results per pool (again fix the souporcell assignment for pools without SNP data)
samples_without_SNP <- c("4456", "4452", "4851", "4389", "4814")

demuxafy_metadata <- list()

for(capture_id in capture_ids){
    result <- read.csv(paste0(path, capture_id, "/combined_results_w_combined_assignments.tsv"), sep="\t", header=TRUE)

    index_threshold <- grep(threshold, colnames(result))
    index_droplet <- grep(droplet, colnames(result))

    samples <- unique(result$Souporcell_Cluster)
    for(sample in samples){
        if(sample%in%samples_without_SNP){
            sample
            find <- grep(sample, result$Souporcell_Cluster)
            result[find, index_threshold] <- sample
            result[find, index_droplet] <- "singlet"
        }
    }
    #write.table(result, file=paste0(path, capture_id, "/combined_results_w_combined_assignments_edited.tsv"), 
            #append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

   
    # write metadata for scRNAseq run (barcode + final assignment) 
    result_sc <- result
    meta <- metadata %>% filter(CAPTURE_ID == paste(capture_id))
    result_sc[,index_threshold] <- meta$SAMPLE_ID[match(result_sc[,index_threshold], meta$CEL_sample)]
    result_sc[,index_threshold] <- ifelse(is.na(result_sc[,index_threshold] ), result_sc[,index_droplet], result_sc[,index_threshold] )
    scRNAseq_metadata <- result_sc %>% select(Barcode, paste(threshold), paste(droplet))
    # change colnames: barcode, assignment (sampleID or doublet), status (singlet, doublet or unassigned)
    colnames(scRNAseq_metadata) <- c("barcode", "assignment", "status")
    scRNAseq_metadata$status[grep("singlet", scRNAseq_metadata$assignment)] <- "unassigned"
    #write.table(scRNAseq_metadata, file=paste0(path, capture_id, "/scRNAseq_workflow_config.tsv"), 
            #append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

    
    result[,index_threshold] <- meta$SAMPLE_ID[match(result[,index_threshold], meta$CEL_sample)]
    result[,index_threshold] <- ifelse(is.na(result[,index_threshold] ), result[,index_droplet], result[,index_threshold] )
    result <- result %>% select(Barcode, paste(threshold))
    result$pool <- paste(capture_id)
    demuxafy_metadata[[capture_id]] <- result
}

# unlist the list of list into a data frame 
df <- as.data.frame(do.call(rbind, demuxafy_metadata))
rownames(df) <- NULL

#used to be /demuxafy_metadata_final.tsv
write.table(df, file=paste0(path, "report/demuxafy_metadata.tsv"), 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

all <- list()
singlets <- list()

# make tables with all cell counts, singlets, unassigned, doublets 
for(capture_id in capture_ids){
  # cell count table 
    table <- result %>% select(paste(threshold))
    table <- as.data.frame(table(table))
  
    sum <- sum(table$Freq)
    table$fraction_all_cells <- paste0(round((table$Freq)/sum * 100), "%")
    colnames(table) <- c("sample", "number_of_cells", "fraction_of_all_cells_captured")
  
    # singlets
    table_singlets <- filter(paste(threshold) != "doublet") %>% filter(paste(threshold) != "singlet") %>% select(paste(threshold))
    table_singlets <- as.data.frame(table(table_singlets))
    sum_singlets <- sum(table_singlets$Freq)
  
    #table_singlets_no_PBMC <- table_singlets %>% filter(MajoritySinglet_Individual_Assignment != "PBMC-Donor-6")
    table_singlets_no_PBMC <- table_singlets %>% filter(table_singlets!="PBMC_spike")
    sum_singlets_no_PBMC <-  sum(table_singlets_no_PBMC$Freq)
    table_singlets_no_PBMC$fraction_singlets_no_PBMC <- paste0(round((table_singlets_no_PBMC$Freq)/sum_singlets_no_PBMC * 100), "%")
  
    #get back to singlet table
    table_singlets$fraction_singlets <- paste0(round((table_singlets$Freq)/sum_singlets * 100), "%")
    #table$fraction_singlets <- c(fraction_singlets, "NA", "NA")
  
    colnames(table_singlets) <- c("sample", "number_of_cells", "fraction_of_all_singlets_captured")
    colnames(table_singlets_no_PBMC) <- c("sample", "number_of_cells", "fraction_of_singlets_captured_without_PBMC")
  
    all_table_singlets <- full_join(table_singlets, table_singlets_no_PBMC, by = c("sample", "number_of_cells"))
    all_table <- full_join(table, all_table_singlets, by = c("sample", "number_of_cells"))
    all_table$pool <- paste(capture_id)
    samples_found <- nrow(table_singlets)
  
    doublets <- table %>% filter(sample == "doublet") %>% pull("number_of_cells")
    unassigned <- table %>% filter(sample == "unassigned") %>% pull("number_of_cells")
    PBMCs <- table %>% filter(sample == "PBMC_spike") %>% pull("number_of_cells") 
  
    all_table_singlets$pool <- paste(capture_id)
    reactions <- all_table_singlets$pool
    reactions <- gsub("_repeat", "", reactions, ignore.case = T)
    reactions <- substr(reactions,1,nchar(reactions)-1) 
    all_table_singlets$reactions <- reactions
    all_table_singlets$reactions <- gsub("Atlas_Pool2$", "Atlas_Pool24", all_table_singlets$reactions)
    all_table_singlets$reactions <- gsub("Atlas_buffy_$", "Atlas_buffy", all_table_singlets$reactions)
  
    reactions <- all_table$pool
    reactions <- gsub("_repeat", "", reactions, ignore.case = T)
    reactions <- substr(reactions,1,nchar(reactions)-1) 
    all_table$reactions <- reactions
    all_table$reactions <- gsub("Atlas_Pool2$", "Atlas_Pool24", all_table$reactions)
    all_table$reactions <- gsub("Atlas_buffy_$", "Atlas_buffy", all_table$reactions)
    all[[capture_id]] <- all_table
    singlets[[capture_id]] <- all_table_singlets
}

# unlist the lists of all cells and only singlets
df_all <- as.data.frame(do.call(rbind, all))
rownames(df_all) <- NULL
df_singlets <- as.data.frame(do.call(rbind, singlets))
rownames(df_singlets) <- NULL

write.table(df_all, file=paste0(path, "report/all_table.tsv"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(df_singlets, file=paste0(path, "report/all_table_singlets.tsv"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)