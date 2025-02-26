source("scripts/ckn_utils/ckn_utils_load_peaks.R")

# all_H3_macs2peaks <- load_h3_macs2peaks()
# all_H3_gopeaks <- load_h3_gopeaks()
# 
# H3K4me3_macs2peaks <- all_H3_macs2peaks[grepl("H3K4me3", names(all_H3_macs2peaks))]
# names(H3K4me3_macs2peaks) <- paste0(names(H3K4me3_macs2peaks), "_macs2")
# H3K27me3_macs2peaks <- all_H3_macs2peaks[grepl("H3K27me3", names(all_H3_macs2peaks))]
# names(H3K27me3_macs2peaks) <- paste0(names(H3K27me3_macs2peaks), "_macs2")
# H3K4me3_gopeaks <- all_H3_gopeaks[grepl("H3K4me3", names(all_H3_gopeaks))]
# H3K27me3_gopeaks <- all_H3_gopeaks[grepl("H3K27me3", names(all_H3_gopeaks))]
# 
# all_H3K4me3_peaks <- c(H3K4me3_macs2peaks, H3K4me3_gopeaks)
# all_H3K27me3_peaks <- c(H3K27me3_macs2peaks, H3K27me3_gopeaks)

#
peaks <- all_RAJI_peaks[[names(all_RAJI_peaks)[1]]]
# seqlevelsStyle(peaks) <- "NCBI"

plotRanges(data = peaks, params = params_i,
           y = "0.25b", height = 0.05,
           linecolor = NA, fill =  "#4A6EA3", collapse = TRUE)
plotText(label = names(all_RAJI_peaks)[1], x = 0.95, y = "-0.025b", fontsize = fontsize_val, just = "right")

for (peakset in names(all_RAJI_peaks)[2]) {
  message("# ", peakset)
  
  peaks <- all_RAJI_peaks[[peakset]]
  # seqlevelsStyle(peaks) <- "NCBI"
  
  plotRanges(data = peaks, params = params_i,
             y = "0.02b", height = 0.05,
             linecolor = NA, fill =  "#4A6EA3", collapse = TRUE)
  plotText(label = peakset, x = 0.95, y = "-0.025b", fontsize = fontsize_val, just = "right")
}

for (peakset in names(all_U2OS_peaks)) {
  message("# ", peakset)
  
  peaks <- all_U2OS_peaks[[peakset]]
  # seqlevelsStyle(peaks) <- "NCBI"
  
  plotRanges(data = peaks, params = params_i,
             y = "0.02b", height = 0.05,
             linecolor = NA, fill =  "#DAA520", collapse = TRUE)
  plotText(label = peakset, x = 0.95, y = "-0.025b", fontsize = fontsize_val, just = "right")
}



# 
# source("scripts/ckn_utils/ckn_utils_H3K27ac_diffbind.R")
# source("scripts/ckn_utils/ckn_utils_cofactor_diffbind.R")
# source("scripts/ckn_utils/ckn_utils_ATAC_diffbind.R")
# 
# #
# plotRanges(data = "/Users/chris/Desktop/20230130_gr_project/output/chip-pipeline-GR-GRCh38_SE/peak_call_withWCE/A549_GR_1h_consensus.GR_peaks.narrowPeak.stdchr.bed", params = params_i,
#            y = "0.25b", height = 0.05,
#            linecolor = NA, fill =  "#303030", collapse = TRUE)
# plotText(label = "GR_1h", x = 0.5, y = "-0.025b", fontsize = fontsize_val-2, just = "right")
# 
# #
# diffbind_cofactor <- load_diffbind_MB_peaks()
# 
# plotRanges(data = diffbind_cofactor[["MBcomb5pT_UP"]], params = params_i,
#            y = "0.05b", height = 0.05,
#            linecolor = NA, fill =  "#4A6EA3", collapse = TRUE)
# plotText(label = "COF_UP", x = 0.5, y = "-0.025b", fontsize = fontsize_val-2, just = "right")
# 
# plotRanges(data = diffbind_cofactor[["MBcomb5pT_UNBIASED"]], params = params_i,
#            y = "0.05b", height = 0.05,
#            linecolor = NA, fill =  "#4A6EA3", collapse = TRUE)
# plotText(label = "COF_UNBIASED", x = 0.5, y = "-0.025b", fontsize = fontsize_val-2, just = "right")
# 
# plotRanges(data = diffbind_cofactor[["MBcomb5pT_DOWN"]], params = params_i,
#            y = "0.05b", height = 0.05,
#            linecolor = NA, fill =  "#4A6EA3", collapse = TRUE)
# plotText(label = "COF_DOWN", x = 0.5, y = "-0.025b", fontsize = fontsize_val-2, just = "right")
# 
# #
# H3K27ac_UP <- load_H3K27ac_UP_overtime_peaks(filter_nb = "filter2", FC = 0.75)
# H3K27ac_UNBIASED <- load_H3K27ac_UNBIASED_overtime_peaks(filter_nb = "filter2", FC = 0.75)
# H3K27ac_DOWN <- load_H3K27ac_DOWN_overtime_peaks(filter_nb = "filter2", FC = 0.75)
# 
# plotRanges(data = H3K27ac_UP, params = params_i,
#            y = "0.05b", height = 0.05,
#            linecolor = NA, fill =  "#00441B", collapse = TRUE)
# plotText(label = "H3K27ac_UP", x = 0.5, y = "-0.025b", fontsize = fontsize_val-2, just = "right")
# 
# plotRanges(data = H3K27ac_UNBIASED, params = params_i,
#            y = "0.05b", height = 0.05,
#            linecolor = NA, fill =  "#00441B", collapse = TRUE)
# plotText(label = "H3K27ac_UNBIASED", x = 0.5, y = "-0.025b", fontsize = fontsize_val-2, just = "right")
# 
# plotRanges(data = H3K27ac_DOWN, params = params_i,
#            y = "0.05b", height = 0.05,
#            linecolor = NA, fill =  "#00441B", collapse = TRUE)
# plotText(label = "H3K27ac_DOWN", x = 0.5, y = "-0.025b", fontsize = fontsize_val-2, just = "right")
# 
# #
# diffbind_ATAC <- load_diffbind_ATAC_peaks(timepoint_list = "1h", filter_nb = "custom3", FC = 0.75)
# 
# plotRanges(data = diffbind_ATAC[["ATAC_1h_UP_FC0p75"]], params = params_i,
#            y = "0.05b", height = 0.05,
#            #linecolor = NA, fill =  "#C0392B", collapse = TRUE)
#            linecolor = NA, fill =  "#4CA189", collapse = TRUE)
# plotText(label = "ATAC_UP", x = 0.5, y = "-0.025b", fontsize = fontsize_val-2, just = "right")
# 
# plotRanges(data = diffbind_ATAC[["ATAC_1h_UNBIASED_FC0p75"]], params = params_i,
#            y = "0.05b", height = 0.05,
#            linecolor = NA, fill =  "#4CA189", collapse = TRUE)
# plotText(label = "ATAC_UNBIASED", x = 0.5, y = "-0.025b", fontsize = fontsize_val-2, just = "right")
# 
# plotRanges(data = diffbind_ATAC[["ATAC_1h_DOWN_FC0p75"]], params = params_i,
#            y = "0.05b", height = 0.05,
#            linecolor = NA, fill =  "#4CA189", collapse = TRUE)
# plotText(label = "ATAC_DOWN", x = 0.5, y = "-0.025b", fontsize = fontsize_val-2, just = "right")