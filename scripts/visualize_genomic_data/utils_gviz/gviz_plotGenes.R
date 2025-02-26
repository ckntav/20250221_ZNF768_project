##### Plot gene track
plotGenes(params = params_i,
          chrom = gsub("chr", "", chr_i),
          y = "0.025b", height = 0.5,
          fontcolor = c("#404040", "#404040"),
          fill = c("#404040", "#404040"),
          stroke = 0.01,
          # bg = "red",
          assembly = hg38_ensembl104,
          fontsize = fontsize_val)

# ##### Plot chromatin states
# #
# chromstate_0h_path <- "input/chromatin_states/A549_ChromHMM_18state_model_0h_ENCFF418WHV.bed.gz"
# chromstate_0h <- rtracklayer::import(chromstate_0h_path) %>% as.data.frame
# 
# #
# # chromstate_1h_path <- "input/chromatin_states/A549_ChromHMM_18state_model_1h_ENCFF052NXZ.bed.gz"
# # chromstate_1h <- rtracklayer::import(chromstate_1h_path) %>% as.data.frame
# 
# #
# chromstate_2h_path <- "input/chromatin_states/A549_ChromHMM_18state_model_2h_ENCFF246IPY.bed.gz"
# chromstate_2h <- rtracklayer::import(chromstate_2h_path) %>% as.data.frame
# 
# #
# chromstate_4h_path <- "input/chromatin_states/A549_ChromHMM_18state_model_4h_ENCFF910RII.bed.gz"
# chromstate_4h <- rtracklayer::import(chromstate_4h_path) %>% as.data.frame
# 
# 
# #
# df_chromstate <- chromstate_0h %>%
#   dplyr::select(name, itemRgb) %>% 
#   distinct()
# 
# color_chromstate <- df_chromstate$itemRgb
# chromstate_0h$name <- factor(chromstate_0h$name, levels = df_chromstate$name)
# # chromstate_1h$name <- factor(chromstate_1h$name, levels = df_chromstate$name)
# chromstate_2h$name <- factor(chromstate_2h$name, levels = df_chromstate$name)
# chromstate_4h$name <- factor(chromstate_4h$name, levels = df_chromstate$name)
# 
# 
# 
# #
# plotRanges(data = chromstate_0h, params = params_i,
#            y = "0.025b", height = 0.05,
#            fill = colorby(column = "name",
#                           palette = colorRampPalette(color_chromstate)), collapse = TRUE)
# plotText(label = "0h", x = 0.45, y = "-0.025b", fontsize = fontsize_val-2, just = "right")
# 
# #
# # plotRanges(data = chromstate_1h, params = params_i,
# #            y = "0.025b", height = 0.05,
# #            fill = colorby(column = "name",
# #                           palette = colorRampPalette(color_chromstate)), collapse = TRUE)
# # plotText(label = "1h", x = 0.45, y = "-0.025b", fontsize = fontsize_val-2, just = "right")
# 
# #
# plotRanges(data = chromstate_2h, params = params_i,
#            y = "0.025b", height = 0.05,
#            fill = colorby(column = "name",
#                           palette = colorRampPalette(color_chromstate)), collapse = TRUE)
# plotText(label = "2h", x = 0.45, y = "-0.025b", fontsize = fontsize_val-2, just = "right")
# 
# #
# plotRanges(data = chromstate_4h, params = params_i,
#            y = "0.025b", height = 0.05,
#            fill = colorby(column = "name",
#                           palette = colorRampPalette(color_chromstate)), collapse = TRUE)
# plotText(label = "4h", x = 0.45, y = "-0.025b", fontsize = fontsize_val-2, just = "right")

##### Plot genome label
plotGenomeLabel(params = params_i,
                # chrom = paste0("chr", params_i$chrom),
                # assembly = ensembl104,
                y = "0.05b", scale = "bp",
                fontcolor = "#404040", linecolor = "#404040",
                fontsize = fontsize_val,
                boxWidth = 0.2)


# chromstate_30m_path <- "input/chromatin_states/A549_ChromHMM_18state_model_30m_ENCFF113TCU.bed.gz"
# chromstate_3h_path <- "input/chromatin_states/A549_ChromHMM_18state_model_3h_ENCFF513UFQ.bed.gz"
# chromstate_5h_path <- "input/chromatin_states/A549_ChromHMM_18state_model_5h_ENCFF662GGJ.bed.gz"
# chromstate_6h_path <- "input/chromatin_states/A549_ChromHMM_18state_model_6h_ENCFF845TIM.bed.gz"
# chromstate_30m <- rtracklayer::import(chromstate_30m_path) %>% as.data.frame
# chromstate_3h <- rtracklayer::import(chromstate_3h_path) %>% as.data.frame
# chromstate_5h <- rtracklayer::import(chromstate_5h_path) %>% as.data.frame
# chromstate_6h <- rtracklayer::import(chromstate_6h_path) %>% as.data.frame
# chromstate_30m$name <- factor(chromstate_30m$name, levels = df_chromstate$name)
# chromstate_3h$name <- factor(chromstate_3h$name, levels = df_chromstate$name)
# chromstate_5h$name <- factor(chromstate_5h$name, levels = df_chromstate$name)
# chromstate_6h$name <- factor(chromstate_6h$name, levels = df_chromstate$name)

# #
# plotRanges(data = chromstate_30m, params = params_i,
#            y = "0.1b", height = 0.15,
#            fill = colorby(column = "name",
#                           palette = colorRampPalette(color_chromstate)), collapse = TRUE)
# plotText(label = "30m", x = 0.45, y = "-0.05b", fontsize = fontsize_val-2, just = "right")
# 
# #
# plotRanges(data = chromstate_3h, params = params_i,
#            y = "0.1b", height = 0.15,
#            fill = colorby(column = "name",
#                           palette = colorRampPalette(color_chromstate)), collapse = TRUE)
# plotText(label = "3h", x = 0.45, y = "-0.05b", fontsize = fontsize_val-2, just = "right")
# 
# 
# #
# plotRanges(data = chromstate_5h, params = params_i,
#            y = "0.1b", height = 0.15,
#            fill = colorby(column = "name",
#                           palette = colorRampPalette(color_chromstate)), collapse = TRUE)
# plotText(label = "5h", x = 0.45, y = "-0.05b", fontsize = fontsize_val-2, just = "right")
# 
# #
# plotRanges(data = chromstate_6h, params = params_i,
#            y = "0.1b", height = 0.15,
#            fill = colorby(column = "name",
#                           palette = colorRampPalette(color_chromstate)), collapse = TRUE)
# plotText(label = "6h", x = 0.45, y = "-0.05b", fontsize = fontsize_val-2, just = "right")