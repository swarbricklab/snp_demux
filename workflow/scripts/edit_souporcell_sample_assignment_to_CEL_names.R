# edit soupocell assignment to extract the CEL based sample names and add them to the souporcell assignments (prior to running Combine_results.R)
dir = "/Users/beatak/external/Cluster230201/"
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
  clusters_file <- read.csv(paste0(dir, "results/BrCaAtlas/", p, "/annotation/vireosnp/clusters_annotated.tsv"), sep="\t", header=TRUE)
  clusters_file_singlet <- clusters_file %>% filter(status == "singlet")
  clusters_file_singlet$assignment <- meta$CEL_sample[match(clusters_file_singlet$assignment, meta$SAMPLE_ID)]
  clusters_file$assignment[match(clusters_file_singlet$barcode, clusters_file$barcode)] <- clusters_file_singlet$assignment
  #clusters_file <- clusters_file %>% replace_na(list(assignment == 'unasigned'))
  write.table(clusters_file, file=paste0(dir, "results/BrCaAtlas/", p, "/annotation/vireosnp/clusters_annotated_to_match_demuxafy.tsv"), 
              append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# next run Combine_results.R