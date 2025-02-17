# ===================
# title: "demux_souporcell.generate_cell_barcode_lists"
# description: extracts cell-barcode lists from souporcell generated sample demux
# project: "workflow-single-cell-rna-seq"
# author: "Daniel Roden"
# email: "d.roden@garvan.org.au / dlroden@gmail.com"
# company: "Garvan Institute"
# created: 22/09/2020
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
result_file_Demuxlet <- snakemake@input[["result_Demuxlet"]]
result_file_Freemuxlet <- snakemake@input[["result_Freemuxlet"]]
result_file_Souporcell <- snakemake@input[["result_Souporcell"]]
result_file_Vireo <- snakemake@input[["result_Vireo"]]

message("Demuxlet file: ", result_file_Demuxlet)
message("Freemuxlet file: ", result_file_Freemuxlet)
message("Souporcell file: ", result_file_Souporcell)
message("Vireo file: ", result_file_Vireo)

message("OUTPUT")
file_upset_plot_pool <- snakemake@output[["upset_plot_pool"]]
message("upset plot for pool: ", file_barcode_list)
file_upset_plot_per_sample <- snakemake@output[["upset_plot_per_sample"]]
message("upset plot per sample in a pool: ", file_cell_counts)


message("WILDCARDS")
capture_id <- snakemake@wildcards[["capture_id"]]
message("capture_id: ", capture_id)

message("----------------------------------")
message("EXTRACT SINGLET CELLs in SAMPLEs") # ----

result_Demuxlet = read.table(result_file_Demuxlet, header=T)
result_Freemuxlet = read.table(result_file_Freemuxlet, header=T)
result_Souporcell = read.csv(result_file_Souporcell, sep = "\t", header = T)
result_Vireo = read.csv(result_file_Vireo, sep = "\t", header = T)

singlets_Demuxlet <- result_Demuxlet[which(result_Demuxlet$DROPLET.TYPE =="SNG"),"BARCODE"]
#add freemuxlet
singlets_Souporcell <- result_Souporcell[which(result_Souporcell$status == "singlet"),"barcode"]
singlets_Vireo <- result_Vireo[which(result_Vireo$donor_id != "unassigned"),]
singlets_Vireo <- singlets_Vireo[which(singlets_Vireo$donor_id != "doublet"),"cell"]

message("EXTRACTED ", nrow(singlets_Demuxlet), " singlets from Demuxlet")
#message("EXTRACTED ", nrow(singlets_Freemuxlet), " singlets from Freemuxlet")
message("EXTRACTED ", nrow(singlets_Souporcell), " singlets from Souporcell")
message("EXTRACTED ", nrow(singlets_Vireo), " singlets from Vireo")

listInput <- list(Demuxlet = singlets_Demuxlet, Freemuxlet = singlets_Freemuxlet, Souporcell = singlets_Souporcell, Vireo = singlets_Vireo)

pdf(file = file_upset_plot_pool, width = 4, height = 4)
upset(fromList(listInput), order.by = "freq")
dev.off()




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
