library(tidyverse)

#
wdir_path <- "/Users/chris/Desktop/20250221_ZNF768_project"
peaks_antibody_dir <- file.path(wdir_path, "output", "GEO_peak_GSE76496")
chromsizeshg38 <- file.path(wdir_path, "input/genome_annot/hg38.chrom.sizes")

#
sample_list <- c("HEK293_ZNF768_rep1.bed")

for (sample in sample_list) {
  message("# ", sample)
  
  input_bed <- sample
  input_path <- file.path(peaks_antibody_dir, input_bed)
  input_path_sorted <- str_replace(input_path, pattern = "\\.bed", replacement = "\\.sorted.bed")
  
  output_bb <- str_replace(input_path_sorted, pattern = "bed", replacement = "bb")
  # message(output_bb)
  # output_path <- file.path(peaks_antibody_dir, output_bb)
  
  # sort 
  call_sort <- paste("sort -k1,1 -k2,2n",
                     input_path, ">",
                     input_path_sorted)
  
  system(call_sort)
  
  call_bedToBigBed <- paste("/Users/chris/Documents/software/bedToBigBed",
                            "-type=bed3+3",
                            input_path_sorted, chromsizeshg38, output_bb)
  # message(call_bedToBigBed)
  system(call_bedToBigBed)
}

# 
# # for (i in (which(df$antibody != "WCE"))) {
# for (i in 1:nrow(df)) {
#   sample_name <- df$sample_name[i]
#   antibody_factor <- df$antibody_factor[i]
#   condition_i <- df$condition[i]
#   treatment_i <- df$treatment[i]
#   rep_i <- df$replicate[i]
#   # message("# ", sample_name)
#   
#   #
#   basename <- paste("LAPC4", condition_i, antibody_factor, treatment_i, rep_i, sep = "_")
#   
#   input_bed <- paste0(basename, ".",  antibody_factor, "_peaks.narrowPeak.stdchr.bed")
#   input_path <- file.path(peaks_antibody_dir, basename,  antibody_factor, input_bed)
#   output_bb <- str_replace(input_bed, pattern = "bed", replacement = "bb")
#   # message(output_bb)
#   output_path <- file.path(tracks_dir, output_bb)
#   
#   # 
#   call_bedToBigBed <- paste("/Users/chris/Documents/software/bedToBigBed",
#                             "-type=bed3+3",
#                             input_path, chromsizeshg38, output_path)
#   # message(call_bedToBigBed)
#   system(call_bedToBigBed)
# }
