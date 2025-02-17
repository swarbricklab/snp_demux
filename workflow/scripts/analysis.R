results from our souporcell, annotated: singlet	5%_PBMC_spike_ins, singlet	4666, etc 
/directflow/TumourProgressionGroupTemp/projects/tmp/beakie/brca_atlas/brca_atlas-data/soup-or-cell/results/BrCaAtlas/Atlas_Pool_9a/annotation/vireosnp/clusters_annotated.tsv 
/directflow/TumourProgressionGroupTemp/projects/tmp/beakie/brca_atlas/brca_atlas-data/results/BrCaAtlas/Atlas_Pool_9a/annotation/vireosnp

which looks like the output from demuxafy-souporcell, but need to change the names to match the SNP based names 
# file called clusters.tsv here: /directflow/TumourProgressionGroupTemp/projects/tmp/beakie/demuxafy/results_Souporcell/clusters.tsv

dir = "/Users/beatak/external/Cluster230131lastlast/"
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

# here the combine_results.R can be run - part of the demuxafy 

# read the combined results per pool (again fix the souporcell assignment for pools without SNP data)
samples_without_SNP <- c("4456", "4452", "4851", "4389", "4814")
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

# write metadata for scRNAseq run (barcode + final assignment) 
for(p in pools){
  result <- read.csv(paste0(dir, "results/demuxafy/combined/BrCaAtlas/", p, "/combined_results_w_combined_assignments_edited.tsv"), sep="\t", header=TRUE)
  meta <- metadata %>% filter(CAPTURE_ID == p)
  result$MajoritySinglet_Individual_Assignment <- meta$SAMPLE_ID[match(result$MajoritySinglet_Individual_Assignment, meta$CEL_sample)]
  result$MajoritySinglet_Individual_Assignment <- ifelse(is.na(result$MajoritySinglet_Individual_Assignment ), result$MajoritySinglet_DropletType, result$MajoritySinglet_Individual_Assignment )
  scRNAseq_metadata <- result %>% select(Barcode, MajoritySinglet_Individual_Assignment)
  write.table(scRNAseq_metadata, file=paste0(dir, "results/demuxafy/combined/BrCaAtlas/", p, "/scRNAseq_workflow_metadata.tsv"), 
              append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# read the metadata file and add pool name 
demuxafy_metadata<- list()
for(p in pools){
  result <- read.csv(paste0(dir, "results/demuxafy/combined/BrCaAtlas/", p, "/combined_results_w_combined_assignments_edited.tsv"), sep="\t", header=TRUE)
  meta <- metadata %>% filter(CAPTURE_ID == p)
  result$MajoritySinglet_Individual_Assignment <- meta$SAMPLE_ID[match(result$MajoritySinglet_Individual_Assignment, meta$CEL_sample)]
  result$MajoritySinglet_Individual_Assignment <- ifelse(is.na(result$MajoritySinglet_Individual_Assignment ), result$MajoritySinglet_DropletType, result$MajoritySinglet_Individual_Assignment )
  result <- result %>% select(Barcode, MajoritySinglet_Individual_Assignment)
  result$pool <- p 
  demuxafy_metadata[[p]] <- result
}

# unlist the list of list into a table 
df <- as.data.frame(do.call(rbind, demuxafy_metadata))
rownames(df) <- NULL
write.table(df, file=paste0(dir, "results/demuxafy/combined/BrCaAtlas/demuxafy_metadata_final.tsv"), 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

all_cells = 0 
all_singlets = 0 
all_table <- data.frame()
all_table_singlets <- data.frame()
all_doublets = 0
all_unassigned = 0 
all_singlets_no_PBMC = 0 
all_PBMCs = 0

all <- list()
singlets <- list()

# make tables with all cell counts, singlets, unassigned, doublets 
for(pool in pools){
  
  # cell count table 
  result <- read.csv(paste0(dir, "results/demuxafy/combined/BrCaAtlas/", pool, "/combined_results_w_combined_assignments_edited.tsv"), sep="\t", header=TRUE)
  meta <- metadata %>% filter(CAPTURE_ID == pool)
  result$MajoritySinglet_Individual_Assignment <- meta$SAMPLE_ID[match(result$MajoritySinglet_Individual_Assignment, meta$CEL_sample)]
  result$MajoritySinglet_Individual_Assignment <- ifelse(is.na(result$MajoritySinglet_Individual_Assignment ), result$MajoritySinglet_DropletType, result$MajoritySinglet_Individual_Assignment )
  table <- result %>% select(MajoritySinglet_Individual_Assignment)
  table <- as.data.frame(table(table))
  
  sum <- sum(table$Freq)
  table$fraction_all_cells <- paste0(round((table$Freq)/sum * 100), "%")
  colnames(table) <- c("sample", "number_of_cells", "fraction_of_all_cells_captured")
  #table <- filter(table, number_of_cells > 3)
  
  # singlets
  table_singlets <- result %>% filter(MajoritySinglet_Individual_Assignment != "doublet") %>% select(MajoritySinglet_Individual_Assignment)
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
  
  #table %>% filter(number_of_cells > 5)
  # this was removing PBMC in some cases, so I decided to keep it
  #table_singlets <- filter(table_singlets, Freq > 3)
  
  all_table_singlets <- full_join(table_singlets, table_singlets_no_PBMC, by = c("sample", "number_of_cells"))
  all_table <- full_join(table, all_table_singlets, by = c("sample", "number_of_cells"))
  all_table$pool <- pool
  samples_found <- nrow(table_singlets)
  
  doublets <- table %>% filter(sample == "doublet") %>% pull("number_of_cells")
  unassigned <- table %>% filter(sample == "unassigned") %>% pull("number_of_cells")
  PBMCs <- table %>% filter(sample == "PBMC_spike") %>% pull("number_of_cells") 
  
  all_table_singlets$pool <- pool
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
  all[[pool]] <- all_table
  singlets[[pool]] <- all_table_singlets
}

# unlist the lists of all cells and only singlets
df_all <- as.data.frame(do.call(rbind, all))
rownames(df_all) <- NULL
df_singlets <- as.data.frame(do.call(rbind, singlets))
rownames(df_singlets) <- NULL

dir_save = "/Users/beatak/external/Cluster230131lastlast/results/demuxafy/report/all/"
#dir.create(dir_save, showWarnings = TRUE, recursive = FALSE)
write.table(df_all, file=paste0(dir_save, "all_table_230131.tsv"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(df_singlets, file=paste0(dir_save, "all_table_singlets_230131.tsv"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  


metadata_subtype <- read.csv("/Users/beatak/Dropbox (Garvan)/SwarbrickLab Team Folder/Beata_Kiedik/BrCaAtlas/metadata/metadata_brca_updated_subtypes.csv")

# for re-capturing, per pool:
df_summary <- as.data.frame(df_all %>%
                              group_by(sample, reactions) %>% 
                              #group_by(sample) %>% 
                              summarise(total_number_of_cells = sum(number_of_cells)))

pool_name_GWCCG <- as.data.frame(unique(metadata_subtype$Pool))
pool_name_GWCCG$pool_name_GWCCG <- c("Atlas_buffy", "Atlas_Pool_11", "Atlas_Pool_CD45_P2", 
                                     "Atlas_CD45_Pool3", "Atlas_CD45_Pool4", "Atlas_Pool_1",
                                     "Atlas_Pool_2", "Atlas_Pool_3", "Atlas_Pool04", "Atlas_Pool_5", 
                                     "Atlas_Pool_6","Atlas_Pool_7","Atlas_Pool_8","Atlas_Pool_9","Atlas_Pool_10", 
                                     "Atlas_Pool_12","Atlas_Pool_13","Atlas_Pool_14","Atlas_Pool15","Atlas_Pool16", 
                                     "Atlas_Pool17","Atlas_Pool18","Atlas_Pool19","Atlas_Pool20","Atlas_Pool21",  
                                     "Atlas_Pool22","Atlas_Pool23","Atlas_Pool24")

colnames(pool_name_GWCCG) <- c("Pool", "Pool_GWCCG")
metadata_subtype <- left_join(metadata_subtype, pool_name_GWCCG, by="Pool")
# prepare small meta with reaction name and subtype
meta_pool_subtype_patient <- metadata_subtype[,c("Pool_GWCCG", "subtype", "patient")]
colnames(meta_pool_subtype_patient) <- c("reactions", "subtype", "patient")

df_summary$patient <- gsub("_.*", "", df_summary$sample)

df_summary <- left_join(df_summary, meta_pool_subtype_patient, by=c("patient", "reactions"))

df_500 <- df_summary %>% filter(total_number_of_cells <= 500) %>% filter(sample != "PBMC_spike")  %>% filter(sample != "unassigned")  %>% filter(sample != "doublet") %>% filter(sample != "singlet")
df_500 <- df_500[!duplicated(df_500),]
write.table(df_500, file=paste0("/Users/beatak/Dropbox (Garvan)/SwarbrickLab Team Folder/Beata_Kiedik/BrCaAtlas/metadata/samples_with_less_than_500_cells.csv"), append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

df_1000 <- df_summary %>% filter(total_number_of_cells <= 1000) %>% filter(sample != "PBMC_spike")  %>% filter(sample != "unassigned") %>% filter(sample != "doublet") %>% filter(sample != "singlet")
df_1000 <- df_1000[!duplicated(df_1000),]
write.table(df_1000, file=paste0("/Users/beatak/Dropbox (Garvan)/SwarbrickLab Team Folder/Beata_Kiedik/BrCaAtlas/metadata/samples_with_less_than_1000_cells.csv"), append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

# here
df_summary <- as.data.frame(df_all %>%
                              filter(sample != "unassigned") %>%
                              filter(sample != "doublet") %>%
                              group_by(reactions, sample) %>% 
                              summarise(total_number_of_cells = sum(number_of_cells)))     

df_summary$patient <- gsub("_.*", "", df_summary$sample)

df_summary <- left_join(df_summary, meta_pool_subtype_patient, by=c("patient", "reactions"))

# plot #cells per subtype 
cells_per_pool <- aggregate(df_summary$total_number_of_cells, by=list(subtype=df_summary$subtype), FUN=sum)

cells_per_pool %>% 
  ggplot(aes(x = subtype, y = x, fill=subtype), show.legend = FALSE) +
  geom_bar(stat="identity", position="dodge", width=0.8, size=3) +
  #facet_wrap(~pool, ncol = 10) + 
  #scale_y_continuous(breaks = round(seq(0, max(cells_per_pool$x), by = 2000),1), sec.axis = sec_axis(~ . * 1, breaks = round(seq(0, max(cells_per_pool$x), by = 2000),1)))+
  #geom_hline(yintercept = c(1000, 3000, 5000, 7500, 10000), linetype='dotted', col = 'red') + 
  #annotate("text", x = "PBMC-Donor-6", y = c(1005, 3005, 5005, 7005, 10005), label = c("1000", "3000,", "5000", "7000", "10000"), vjust = -0.5, size = 3) +
  #geom_text(aes(label=number_of_cells), position=position_dodge(width=0.9), vjust=-0.25) + 
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 6)) + 
  ylab("Number of cells")
  #xlab("Sample")

# plot #cells with phyllodes in one group
df_summary$subtype_main <- df_summary$subtype
df_summary$subtype_main[grep("Phyllodes", df_summary$subtype)] <- "Phyllodes"
df_summary$subtype_main[grep("Fibroadenoma", df_summary$subtype)] <- "Phyllodes"
df_summary$subtype_main[grep("Phyllodes", df_summary$subtype)] <- "Phyllodes"

# divide into 4 groups, remove PBMC and Buffy later in the processing
df_summary$subtype_4 <- df_summary$subtype_main
df_summary$subtype_4[grep("Luminal", df_summary$subtype)] <- "Luminal"

df_summary <- df_summary %>% filter(sample != "singlet")

cells_per_pool <- aggregate(df_summary$total_number_of_cells, by=list(subtype=df_summary$subtype_main), FUN=sum)

cells_per_pool %>% 
  ggplot(aes(x = subtype, y = x, fill=subtype), show.legend = FALSE) +
  geom_bar(stat="identity", position="dodge", width=0.8, size=3) +
  #facet_wrap(~pool, ncol = 10) + 
  #scale_y_continuous(breaks = round(seq(0, max(cells_per_pool$x), by = 2000),1), sec.axis = sec_axis(~ . * 1, breaks = round(seq(0, max(cells_per_pool$x), by = 2000),1)))+
  #geom_hline(yintercept = c(1000, 3000, 5000, 7500, 10000), linetype='dotted', col = 'red') + 
  #annotate("text", x = "PBMC-Donor-6", y = c(1005, 3005, 5005, 7005, 10005), label = c("1000", "3000,", "5000", "7000", "10000"), vjust = -0.5, size = 3) +
  #geom_text(aes(label=number_of_cells), position=position_dodge(width=0.9), vjust=-0.25) + 
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 6)) + 
  ylab("Number of cells")
#xlab("Sample")

# number of samples
subtype_table <- as.data.frame(t(table(df_summary$subtype)))
colnames(subtype_table) <- c("Var1", "subtype", "Freq")

subtype_table %>% 
  ggplot(aes(x = subtype, y = Freq, fill=subtype), show.legend = FALSE) +
  geom_bar(stat="identity", position="dodge", width=0.8, size=3) +
  #facet_wrap(~pool, ncol = 10) + 
  #scale_y_continuous(breaks = round(seq(0, max(cells_per_pool$x), by = 2000),1), sec.axis = sec_axis(~ . * 1, breaks = round(seq(0, max(cells_per_pool$x), by = 2000),1)))+
  #geom_hline(yintercept = c(1000, 3000, 5000, 7500, 10000), linetype='dotted', col = 'red') + 
  #annotate("text", x = "PBMC-Donor-6", y = c(1005, 3005, 5005, 7005, 10005), label = c("1000", "3000,", "5000", "7000", "10000"), vjust = -0.5, size = 3) +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) + 
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 6)) + 
  ylab("Number of samples")
#xlab("Sample")

# pie chart
subtype_table <- subtype_table %>%
  arrange(desc(subtype)) %>%
  mutate(lab.ypos = cumsum(Freq) - 0.6*Freq)

ggplot(subtype_table, aes(x = "", y = Freq, fill = subtype)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = Freq), color = "black")+
  #scale_fill_manual(values = mycols) +
  theme_void()

# main subtype
subtype_table <- as.data.frame(t(table(df_summary$subtype_main)))
colnames(subtype_table) <- c("Var1", "subtype", "Freq")

subtype_table %>% 
  ggplot(aes(x = subtype, y = Freq, fill=subtype), show.legend = FALSE) +
  geom_bar(stat="identity", position="dodge", width=0.8, size=3) +
  #facet_wrap(~pool, ncol = 10) + 
  #scale_y_continuous(breaks = round(seq(0, max(cells_per_pool$x), by = 2000),1), sec.axis = sec_axis(~ . * 1, breaks = round(seq(0, max(cells_per_pool$x), by = 2000),1)))+
  #geom_hline(yintercept = c(1000, 3000, 5000, 7500, 10000), linetype='dotted', col = 'red') + 
  #annotate("text", x = "PBMC-Donor-6", y = c(1005, 3005, 5005, 7005, 10005), label = c("1000", "3000,", "5000", "7000", "10000"), vjust = -0.5, size = 3) +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) + 
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 6)) + 
  ylab("Number of samples")
#xlab("Sample")

# pie chart
subtype_table <- subtype_table %>%
  arrange(desc(subtype)) %>%
  mutate(lab.ypos = cumsum(Freq) - 0.6*Freq)

ggplot(subtype_table, aes(x = "", y = Freq, fill = subtype)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = Freq), color = "black")+
  #scale_fill_manual(values = mycols) +
  theme_void()

# 4 subtypes 
df <- df_summary %>% filter(subtype != "PBMC") %>% filter(subtype != "Buffy")
subtype_table <- as.data.frame(t(table(df$subtype_4)))
colnames(subtype_table) <- c("Var1", "subtype", "Freq")

subtype_table <- subtype_table %>%
  arrange(desc(subtype)) %>%
  mutate(lab.ypos = cumsum(Freq) - 0.6*Freq)

ggplot(subtype_table, aes(x = "", y = Freq, fill = subtype)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = Freq), color = "black")+
  #scale_fill_manual(values = mycols) +
  theme_void()

cells_per_pool <- table(x$subtype)
cells_per_pool <- aggregate(df_summary$total_number_of_cells, by=list(subtype=df_summary$subtype, reaction=df_summary$reactions), FUN=sum)

cells_per_pool %>% 
  ggplot(aes(x = reaction, y = x, fill=subtype), show.legend = FALSE) +
  geom_bar(stat="identity", position="dodge", width=0.8, size=3) +
  #facet_wrap(~pool, ncol = 10) + 
  #scale_y_continuous(breaks = round(seq(0, max(cells_per_pool$x), by = 2000),1), sec.axis = sec_axis(~ . * 1, breaks = round(seq(0, max(cells_per_pool$x), by = 2000),1)))+
  #geom_hline(yintercept = c(1000, 3000, 5000, 7500, 10000), linetype='dotted', col = 'red') + 
  #annotate("text", x = "PBMC-Donor-6", y = c(1005, 3005, 5005, 7005, 10005), label = c("1000", "3000,", "5000", "7000", "10000"), vjust = -0.5, size = 3) +
  #geom_text(aes(label=number_of_cells), position=position_dodge(width=0.9), vjust=-0.25) + 
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 6)) + 
  ylab("Number of cells in pool")
#xlab("Sample")

