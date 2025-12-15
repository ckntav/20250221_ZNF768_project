# plotGenomeLabel(params = params_i,
#                 # assembly = ensembl104,
#                 y = 0.3, scale = "bp",
#                 fontcolor = "#404040", linecolor = "#404040",
#                 fontsize = fontsize_val,
#                 boxWidth = 0.2)
# 
# #
# plotSignal(GR_bw[[1]], params = params_i,
#            y = 0.35, height = height_track, scale = TRUE,
#            ymax = GR_scalefactor[[1]],
#            linecolor = "#696969", fill =  "#696969",
#            label = "GR 0h", fontsize = fontsize_val)

RAJI_tracks <- plotSignal(RAJI_bw, params = params_i,
           y = 0.3, height = height_track, scale = TRUE,
           range = c(0, RAJI_maxRange),
           # ymax = scalefactor_i,
           linecolor = "#4A6EA3", fill = "#4A6EA3",
           # label = "RAJI",
           fontsize = fontsize_val)

# annoYaxis(plot = RAJI_tracks, at = c(0, RAJI_maxRange),
#           axisLine = TRUE, fontsize = 5)

plotText(
  label = "RAJI ZNF768 (pooled)", fontsize = 3, fontcolor = "#4A6EA3",
  x = width_val - 0.5, y = 0.3, just = c("right", "top"),
  default.units = "inches")

# RAJI_peaks <- all_RAJI_peaks[["H3K4me3_narrow_st400sl50"]]
# seqlevelsStyle(H3K4me3_peaks) <- "NCBI"
# plotRanges(data = H3K4me3_peaks, params = params_i,
#              y = "0.02b", height = 0.05,
#              linecolor = NA, fill =  "#4A6EA3", collapse = TRUE)
# plotText(label = "H3K4me3_peaks", x = 0.95, y = "-0.025b", fontsize = fontsize_val-2, just = "right")

plotSignal(U2OS_bw, params = params_i,
           y = "0.1b", height = height_track, scale = TRUE,
           range = c(0, U2OS_maxRange),
           # ymax = scalefactor_i,
           linecolor = "#DAA520", fill = "#DAA520",
           # label = "U2OS",
           fontsize = fontsize_val)

plotText(
  label = "U2OS ZNF768 (pooled)", fontsize = 3, fontcolor = "#DAA520",
  x = width_val - 0.5, y = 0.3 + height_track + 0.1, just = c("right", "top"),
  default.units = "inches")

plotSignal(HEK293_bw, params = params_i,
           y = "0.1b", height = height_track, scale = TRUE,
           range = c(0, HEK293_maxRange),
           # ymax = scalefactor_i,
           linecolor = "#C85A54", fill = "#C85A54",
           # label = "U2OS",
           fontsize = fontsize_val)

plotText(
  label = "HEK293 ZNF768 (rep1)", fontsize = 3, fontcolor = "#C85A54",
  x = width_val - 0.5, y = 0.3 + height_track*2 + 0.1*2, just = c("right", "top"),
  default.units = "inches")

# H3K27me3_peaks <- all_H3K27me3_peaks[["H3K27me3_broad_st1000sl200"]]
# seqlevelsStyle(H3K27me3_peaks) <- "NCBI"
# plotRanges(data = H3K27me3_peaks, params = params_i,
#              y = "0.02b", height = 0.05,
#              linecolor = NA, fill =  "#DAA520", collapse = TRUE)
# plotText(label = "H3K27me3_peaks", x = 0.95, y = "-0.025b", fontsize = fontsize_val-2, just = "right")

# Plot RNA 
# for (RNA_sample_name in names(RNA_bw)) {
#   bw_i <- RNA_bw[[RNA_sample_name]]
#   scalefactor_i <- RNA_scalefactor[[RNA_sample_name]]
#   label_i <- gsub("_", " ", RNA_sample_name)
# 
#   plotSignal(bw_i, params = params_i,
#              y = "0.05b", height = height_track, scale = TRUE,
#              ymax = scalefactor_i,
#              linecolor = RNA_colorpalette[[RNA_sample_name]], fill = RNA_colorpalette[[RNA_sample_name]],
#              label = label_i, fontsize = fontsize_val)
# 
#   # readline(prompt="Press [enter] to continue")
# }

# 
# #
# yGR <- list("yes" = "0.05b", "no" = 0.3)
# 
# #
# plotSignal(GR_bw[[1]], params = params_i,
#            y = yGR[[isTAD]], height = height_track, scale = TRUE,
#            ymax = GR_scalefactor[[1]],
#            linecolor = "#696969", fill =  "#696969",
#            label = "GR 0h", fontsize = fontsize_val)
# 
# # Plot GR ChIP-seq
# for (GR_sample_name in names(GR_bw)[2]) {
#   bw_i <- GR_bw[[GR_sample_name]]
#   scalefactor_i <- GR_scalefactor[[GR_sample_name]]
#   label_i <- gsub("_", " ", GR_sample_name)
#   
#   plotSignal(bw_i, params = params_i,
#              y = "0.05b", height = height_track, scale = TRUE,
#              ymax = scalefactor_i,
#              linecolor = "#696969", fill =  "#696969",
#              label = label_i, fontsize = fontsize_val)
#   
#   # readline(prompt="Press [enter] to continue")
# }
# 
# # Plot MED1 ChIP-seq
# for (MED1_sample_name in names(MED1_bw)) {
#   bw_i <- MED1_bw[[MED1_sample_name]]
#   scalefactor_i <- MED1_scalefactor[[MED1_sample_name]]
#   label_i <- gsub("_", " ", MED1_sample_name)
#   
#   plotSignal(bw_i, params = params_i,
#              y = "0.05b", height = height_track, scale = TRUE,
#              ymax = scalefactor_i,
#              linecolor = "#4A6EA3", fill = "#4A6EA3",
#              label = label_i, fontsize = fontsize_val)
#   
#   # readline(prompt="Press [enter] to continue")
# }
# 
# # Plot BRD4 ChIP-seq
# for (BRD4_sample_name in names(BRD4_bw)) {
#   bw_i <- BRD4_bw[[BRD4_sample_name]]
#   scalefactor_i <- BRD4_scalefactor[[BRD4_sample_name]]
#   label_i <- gsub("_", " ", BRD4_sample_name)
#   
#   plotSignal(bw_i, params = params_i,
#              y = "0.05b", height = height_track, scale = TRUE,
#              ymax = scalefactor_i,
#              linecolor = "#E9B65D", fill =  "#E9B65D",
#              label = label_i, fontsize = fontsize_val)
#   
#   # readline(prompt="Press [enter] to continue")
# }
# 
# #
# for (H3K27ac_sample_name in names(H3K27ac_bw)) {
#   bw_i <- H3K27ac_bw[[H3K27ac_sample_name]]
#   scalefactor_i <- H3K27ac_scalefactor[[H3K27ac_sample_name]]
#   label_i <- gsub("_", " ", H3K27ac_sample_name)
#   
#   plotSignal(bw_i, params = params_i,
#              y = "0.05b", height = height_track, scale = TRUE,
#              ymax = scalefactor_i,
#              linecolor = H3K27ac_colorpalette[[H3K27ac_sample_name]], fill = H3K27ac_colorpalette[[H3K27ac_sample_name]],
#              label = label_i, fontsize = fontsize_val)
#   
#   # readline(prompt="Press [enter] to continue")
# }
# 
# # Plot ATAC-seq
# for (ATAC_sample_name in names(ATAC_bw)) {
#   bw_i <- ATAC_bw[[ATAC_sample_name]]
#   scalefactor_i <- ATAC_scalefactor[[ATAC_sample_name]]
#   label_i <- gsub("_", " ", ATAC_sample_name)
#   
#   plotSignal(bw_i, params = params_i,
#              y = "0.05b", height = height_track, scale = TRUE,
#              ymax = scalefactor_i,
#              # linecolor = "#C0392B", fill = "#C0392B",
#              linecolor = "#4CA189", fill = "#4CA189",
#              label = label_i, fontsize = fontsize_val)
  
  # readline(prompt="Press [enter] to continue")
# }