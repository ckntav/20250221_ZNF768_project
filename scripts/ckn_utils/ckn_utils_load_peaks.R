library(GenomicRanges)
library(tidyverse)

#
load_ZNF768_peaks <- function() {
  sample_list <- c("GSM3043267_Raji.R1.bed", "GSM3043268_Raji.R2.bed",
                   "GSM3043270_USOS.R1.bed", "GSM3043271_USOS.R2.bed")
  wdir_path <- "/Users/chris/Desktop/20250221_ZNF768_project"
  peaks_antibody_dir <- file.path(wdir_path, "output", "GEO_peak_GSE111879")
  
  all_peaks <- GRangesList()

  for (sample_i in sample_list) {
    message("# ", sample_i)
    peaks_path <- file.path(peaks_antibody_dir, sample_i)
    peaks <- read_tsv(peaks_path, col_names = FALSE, show_col_types = FALSE)[, 1:3] %>% 
      set_names("seqnames", "start", "end") %>% 
      GRanges()
    # message(" > ", peaks_path)
    message("   > Number of regions : ", length(peaks))
    
    # gather everything
    all_peaks <- append(all_peaks, GRangesList(peaks))
  }
  
  names(all_peaks) <- c("RAJI_ZNF768_rep1", "RAJI_ZNF768_rep2",
                       "U2OS_ZNF768_rep1", "U2OS_ZNF768_rep2")
  message("#####################################")
  message("Available set of regions: ")
  print(names(all_peaks))
  return(all_peaks)
}

# 
# #
# load_LAPC4_peaks <- function(antibodies = c("AR", "ER"),
#                             conditions = c("CTRL", "ovERA"),
#                             treatments = c("V", "T", "TxD", "D"),
#                             reps = c("rep1", "rep2")) {
#   project_dir <- "/Users/chris/Desktop/20240812_LAPC4_ovERa_project"
#   peaks_dir <- "output/chip-pipeline_LAPC4_ovERa_allReps-GRCh38_PE/peak_call_withWCE"
#   
#   all_peaks <- GRangesList()
#   names_all_peaks <- c()
#   
#   for (condition in conditions) {
#     for (antibody in antibodies) {
#       for (treatment in treatments) {
#         for (rep in reps) {
#           sample_name <- paste(sep = "_", "LAPC4", condition, antibody, treatment, rep)
#           message("# ", sample_name)
#           
#           peaks_path <- file.path(project_dir, peaks_dir, sample_name, antibody, paste0(sample_name, ".", antibody, "_peaks.narrowPeak.stdchr.bed"))
#           peaks <- rtracklayer::import(peaks_path)
#           # message(" > ", peaks_path)
#           message("   > Number of regions : ", length(peaks))
#           
#           # gather everything
#           all_peaks <- append(all_peaks, GRangesList(peaks))
#           names_all_peaks <- c(names_all_peaks, sample_name)
#         }
#       }
#     }
#   }
#   names(all_peaks) <- names_all_peaks
#   message("#####################################")
#   message("Available set of regions: ")
#   print(names(all_peaks))
#   return(all_peaks)
# }