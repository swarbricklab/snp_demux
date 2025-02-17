# this script works in R, with cluster mapped to local 
# change souporcell sample ID assignments to match the demuxafy names, which are based on the CEL files

# TODO: fix hardcoded paths
dir = "/Users/beatak/external/Cluster230127/"
#/share/ScratchGeneral/beakie/projects/220802_BrCa_Atlas/brca_atlas-data/soup-or-cell/results/BrCaAtlas
metadata <- read.table("/Users/beatak/Dropbox (Garvan)/SwarbrickLab Team Folder/Beata_Kiedik/BrCaAtlas/metadata/final_souporcell_table_220823.txt", header = TRUE)

# edit metadata to extract the CEL based sample names, 

CEL <- metadata[,7]
# some CEL files dont have the same name pattern = between UKB and Results are is no year, so I added it to be able to extract the sample name in the same way for all CEL files
CEL <- gsub("UKB_RESULTS", "UKB_YYYY_RESULTS", CEL)
CEL <- sub('.*^([^_]+_[^_]+_[^_]+_[^_]+_[^_]+_[^_]+_[^_]+_)', '\\2', CEL)
CEL <- sub("\\..*", "", CEL)
CEL <- gsub("5%_PBMC_spike_ins", "PBMC-Donor-6", CEL)
CEL <- stringr::str_replace_all(CEL, "_", "-")
x <- is.na(CEL) 
CEL[x] <- metadata[x,"SAMPLE_ID"]

metadata$CEL_sample <- CEL
metadata$SAMPLE_ID <- gsub("5%_PBMC_spike_ins", "PBMC_spike", metadata$SAMPLE_ID)  

pools <- unique(metadata$CAPTURE_ID)

for(p in pools){
  meta <- metadata %>% filter(CAPTURE_ID == p)
  # TODO: fix BrCaAtlas project specific paths
  clusters_file <- read.csv(paste0(dir, "results/BrCaAtlas/", p, "/annotation/vireosnp/clusters_annotated.tsv"), sep="\t", header=TRUE)
  clusters_file_singlet <- clusters_file %>% filter(status == "singlet")
  clusters_file_singlet$assignment <- meta$CEL_sample[match(clusters_file_singlet$assignment, meta$SAMPLE_ID)]
  clusters_file$assignment[match(clusters_file_singlet$barcode, clusters_file$barcode)] <- clusters_file_singlet$assignment
  #clusters_file <- clusters_file %>% replace_na(list(assignment == 'unasigned'))
  write.table(clusters_file, file=paste0(dir, "results/BrCaAtlas/", p, "/annotation/vireosnp/clusters_annotated_to_match_demuxafy.tsv"), 
              append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# after combine_results change assignemnts for cells that did not have SNP data 
# these cells were only called by souporcell - need to trust souporcell here :D 
# TO DO fix this crazy for loop to something that makes more sense and is much faster

# TODO: fix hardcoded sampleIDs
samples_without_SNP <- c("4456", "4452", "4851", "4389")
for(p in pools){
  result <- read.csv(paste0(dir, "results/demuxafy/combined/BrCaAtlas/", p, "/combined_results_w_combined_assignments.tsv"), sep="\t", header=TRUE)
  samples <- unique(result$Souporcell_Cluster)
  for(sample in samples){
    if(sample%in%samples_without_SNP){
      sample
      find <- grep(sample, result$Souporcell_Cluster)
      result$MajoritySinglet_Individual_Assignment[find] <- sample
      result$MajoritySinglet_DropletType[find] <- "singlet"
      head(result)
    }
  }
  unique(result$MajoritySinglet_Individual_Assignment)
  write.table(result, file=paste0(dir, "results/demuxafy/combined/BrCaAtlas/", p, "/combined_results_w_combined_assignments_edited.tsv"), 
              append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# generate config/metadata file for scRNAseq workflow
# barcode + assignment
# TO DO  check whether the colnames match what is required in the workflow

demuxafy_metadata<- list()

for(p in pools){
  result <- read.csv(paste0(dir, "results/demuxafy/combined/BrCaAtlas/", p, "/combined_results_w_combined_assignments_edited.tsv"), sep="\t", header=TRUE)
  scRNAseq_metadata <- result %>% filter(MajoritySinglet_DropletType == "singlet") %>% select(Barcode, MajoritySinglet_Individual_Assignment)
  write.table(scRNAseq_metadata, file=paste0(dir, "results/demuxafy/combined/BrCaAtlas/", p, "/scRNAseq_workflow_metadata.tsv"), 
              append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

for(p in pools){
  result <- read.csv(paste0(dir, "results/demuxafy/combined/BrCaAtlas/", p, "/combined_results_w_combined_assignments_edited.tsv"), sep="\t", header=TRUE)
  result <- result %>% select(Barcode, MajoritySinglet_Individual_Assignment)
  result$pool <- p 
  demuxafy_metadata[[p]] <- result
}

demuxafy_metadata_final <- unlist(demuxafy_metadata)
write.table(demuxafy_metadata_final, file=paste0(dir, "results/demuxafy/combined/BrCaAtlas/demuxafy_metadata_final.tsv"), 
              append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)