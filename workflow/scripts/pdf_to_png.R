# ===================
# title: "pdf_to_png"
# description: converts plots in pdf format to png format
# project: "workflow-single-cell-rna-seq - demuxafy"
# author: "Beata Kiedik"
# email: "b.kiedik@garvan.org.au / bkiedik@gmail.com"
# company: "Garvan Institute"
# created: 07/02/2023
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
Singlets_upset_pdf <- snakemake@input[["Singlets_upset_pdf"]]
message("pdf file: ", Singlets_upset_pdf)

message("OUTPUT")
Singlets_upset_png <- snakemake@output[["Singlets_upset_png"]]
message("png file: ", Singlets_upset_png)

message("WILDCARDS")
capture_id <- snakemake@wildcards[["capture_id"]]
message("capture_id: ", capture_id)

message("----------------------------------")
message("Render pdf file to png") # ----

library(png)
library(pdftools)

plot_pdf <- Singlets_upset_pdf
bitmap <- pdf_render_page(plot_pdf, page = 1, dpi = 300)
plot_png <- Singlets_upset_png
png::writePNG(bitmap, plot_png)

message("Converted pdf to png")

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