setwd("/Users/chris/Desktop/20250221_ZNF768_project")

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)  # needed for liftOver

# Read in your peaks
peak_path <- "output/GEO_peak_GSE76496/GSM2026871_ZNF768_rep1.macs.out.txt"
peak_df <- read_tsv(peak_path, skip = 17) %>% 
  dplyr::select(chr, start, end)

# Create GRanges object from hg19 peaks
HEK293_ZNF768_hg19 <- GRanges(peak_df)
HEK293_ZNF768_hg19

# Download the liftOver chain file (only need to do this once)
# chain_file <- "hg19ToHg38.over.chain.gz"
# if (!file.exists(chain_file)) {
#   download.file(
#     url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz",
#     destfile = chain_file
#   )
# }

# Import the chain file
chain_file <- "input/hg19ToHg38.over.chain"
chain <- import.chain(chain_file)

# Perform liftOver
HEK293_ZNF768_hg38 <- liftOver(HEK293_ZNF768_hg19, chain)

# liftOver returns a GRangesList, convert to GRanges
HEK293_ZNF768_hg38 <- unlist(HEK293_ZNF768_hg38)

# Check results
if (TRUE) {
  cat("Original hg19 peaks:", length(HEK293_ZNF768_hg19), "\n")
  cat("Lifted hg38 peaks:", length(HEK293_ZNF768_hg38), "\n")
  cat("Peaks lost in liftOver:", length(HEK293_ZNF768_hg19) - length(HEK293_ZNF768_hg38), "\n")
}

# View the results
HEK293_ZNF768_hg38

# Optional: Export to BED file
export.bed(HEK293_ZNF768_hg38, "output/GEO_peak_GSE76496/HEK293_ZNF768_rep1.bed")



# wdir_path <- "/Users/chris/Desktop/20250221_ZNF768_project"
# peaks_antibody_dir <- file.path(wdir_path, "output", "GEO_peak_GSE111879")
# chromsizeshg38 <- file.path(wdir_path, "input/genome_annot/hg38.chrom.sizes")
# 
# #
# sample_list <- c("GSM3043267_Raji.R1.bed", "GSM3043268_Raji.R2.bed",
#                  "GSM3043270_USOS.R1.bed", "GSM3043271_USOS.R2.bed")
# 
# for (sample in sample_list) {
#   message("# ", sample)
#   
#   input_bed <- sample
#   input_path <- file.path(peaks_antibody_dir, input_bed)
#   input_path_sorted <- str_replace(input_path, pattern = "\\.bed", replacement = "\\.sorted.bed")
#   
#   output_bb <- str_replace(input_path_sorted, pattern = "bed", replacement = "bb")
#   # message(output_bb)
#   # output_path <- file.path(peaks_antibody_dir, output_bb)
#   
#   # sort 
#   call_sort <- paste("sort -k1,1 -k2,2n",
#                      input_path, ">",
#                      input_path_sorted)
#   
#   system(call_sort)
#   
#   call_bedToBigBed <- paste("/Users/chris/Documents/software/bedToBigBed",
#                             "-type=bed3+3",
#                             input_path_sorted, chromsizeshg38, output_bb)
#   # message(call_bedToBigBed)
#   system(call_bedToBigBed)
# }
# # 
# # # for (i in (which(df$antibody != "WCE"))) {
# # for (i in 1:nrow(df)) {
# #   sample_name <- df$sample_name[i]
# #   antibody_factor <- df$antibody_factor[i]
# #   condition_i <- df$condition[i]
# #   treatment_i <- df$treatment[i]
# #   rep_i <- df$replicate[i]
# #   # message("# ", sample_name)
# #   
# #   #
# #   basename <- paste("LAPC4", condition_i, antibody_factor, treatment_i, rep_i, sep = "_")
# #   
# #   input_bed <- paste0(basename, ".",  antibody_factor, "_peaks.narrowPeak.stdchr.bed")
# #   input_path <- file.path(peaks_antibody_dir, basename,  antibody_factor, input_bed)
# #   output_bb <- str_replace(input_bed, pattern = "bed", replacement = "bb")
# #   # message(output_bb)
# #   output_path <- file.path(tracks_dir, output_bb)
# #   
# #   # 
# #   call_bedToBigBed <- paste("/Users/chris/Documents/software/bedToBigBed",
# #                             "-type=bed3+3",
# #                             input_path, chromsizeshg38, output_path)
# #   # message(call_bedToBigBed)
# #   system(call_bedToBigBed)
# # }
