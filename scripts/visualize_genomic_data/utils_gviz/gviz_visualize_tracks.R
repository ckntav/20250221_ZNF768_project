library(GenomicFeatures)

##### Genome assembly
ensembl104_TxDb <- loadDb(file = "input/ensembl/txdb.ensembl104.for.plotgardener.sqlite")
all_transcripts <- transcripts(ensembl104_TxDb, columns = c("tx_name", "gene_id"))
filtered_transcripts <- all_transcripts[all_transcripts$gene_id %in% ensembl_gene_list %>% unlist, ]
mcols(filtered_transcripts)$type <- "transcript"
filtered_txdb <- makeTxDbFromGRanges(filtered_transcripts)

# all_exons <- exons(ensembl104_TxDb, columns = c("tx_name", "gene_id"))
# filtered_exons <- all_exons[all_exons$gene_id %in% ensembl_gene_list %>% unlist, ]
# mcols(filtered_exons)$type <- "exon"
# filtered_txdb <- makeTxDbFromGRanges(filtered_exons)

hg38_ensembl104 <- assembly(
  Genome = "Ensembl.GRCh38.104",
  # TxDb = ensembl104_TxDb,
  TxDb = filtered_txdb,
  OrgDb = "org.Hs.eg.db",
  BSgenome = "BSgenome.Hsapiens.NCBI.GRCh38",
  gene.id.column = "ENSEMBL", display.column = "SYMBOL"
)

##### RAJI
# read RAJI bigwig
readBigwig_RAJI <- function(custom_params) {
  chrom <- custom_params$chrom
  chromstart <- custom_params$chromstart
  chromend <- custom_params$chromend
  message("# Reading RAJI signals at ", paste0("chr", chrom), ":", chromstart, "-", chromend)
  
  #
  bw_dir <- "/Users/chris/Desktop/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/tracks_pooled"
  RAJI_signal <- list()
  bw_filename <- "RAJI_ZNF768_pooled_RPKM.bw"
  bw_filepath <- file.path(bw_dir, bw_filename)
  bw_signal <- readBigwig(file = bw_filepath, params = custom_params)
  RAJI_signal <- bw_signal
}

##### U2OS
# read U2OS bigwig
readBigwig_U2OS <- function(custom_params) {
  chrom <- custom_params$chrom
  chromstart <- custom_params$chromstart
  chromend <- custom_params$chromend
  message("# Reading U2OS signals at ", paste0("chr", chrom), ":", chromstart, "-", chromend)
  
  #
  bw_dir <- "/Users/chris/Desktop/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/tracks_pooled"
  U2OS_signal <- list()
  bw_filename <- "U2OS_ZNF768_pooled_RPKM.bw"
  bw_filepath <- file.path(bw_dir, bw_filename)
  bw_signal <- readBigwig(file = bw_filepath, params = custom_params)
  U2OS_signal <- bw_signal
}

#####
set_maxRange <- function(bw_signal) {
  max_score <- max(bw_signal$score, na.rm = TRUE)
  rounded_max_score <- ceiling(max_score)
  return(rounded_max_score)
}

#####
compute_scalefactor <- function(bw_signal) {
  # Compute the max score for each df
  max_scores <- lapply(bw_signal, function(df) max(df$score, na.rm = TRUE))
  
  # Identify the max score among the max scores
  overall_max_score <- max(unlist(max_scores), na.rm = TRUE)
  
  # Round the overall max score to the ceiling value
  rounded_max_score <- ceiling(overall_max_score/10)*10
  
  # Calculate scale factors for each DataFrame
  scale_factors <- lapply(max_scores, function(score) rounded_max_score/score)
  
  return(scale_factors)
}

##### RNA
# read RNA bigwig
# readBigwig_RNA <- function(custom_params) {
#   chrom <- custom_params$chrom
#   chromstart <- custom_params$chromstart
#   chromend <- custom_params$chromend
#   message("# Reading RNA signals at ", chrom, ":", chromstart, "-", chromend)
# 
#   #
#   bw_dir <- "/Users/chris/Desktop/20240405_bmMSC_HXKX_CnT/output/rnaseq_kallisto"
#   RNA_signal <- list()
#   RNA_samples <- c("bmMSC_RNA_WT_pooled",
#                    "bmMSC_RNA_ME_pooled",
#                    "bmMSC_RNA_JS_pooled")
#   for (sample in RNA_samples) {
#     sample_i <- gsub("_RNA", "", sample)
#     sample_name <- paste(sep = "_", "RNA", sample_i)
#     message("   > ", sample)
#     bw_filename <- paste(sep = "_", sample, "smooth_RPKM.bw")
#     bw_filepath <- file.path(bw_dir, sample, bw_filename)
#     bw_signal <- readBigwig(file = bw_filepath, params = custom_params)
#     RNA_signal[[sample_name]] <- bw_signal
#   }
#   return(RNA_signal)
# }

# determine RNA  color palette
# RNA_colorpalette <- c("#84a59d", "#84a59d", "#4a4e69", "#4a4e69", "#7f4f24", "#7f4f24")
# names(RNA_colorpalette) <- c("RNA_bmMSCWTset1", "RNA_bmMSCWTset2",
#                              "RNA_bmMSCMEset1", "RNA_bmMSCMEset2",
#                              "RNA_bmMSCJSset1", "RNA_bmMSCJSset2")

# RNA_colorpalette <- c("#84a59d", "#4a4e69", "#7f4f24")
# names(RNA_colorpalette) <- c("RNA_bmMSC_WT_pooled", "RNA_bmMSC_ME_pooled", "RNA_bmMSC_JS_pooled")

# 
# ##### H3K27ac
# # read H3K27ac bigwig
# readBigwig_H3K27ac <- function(custom_params) {
#   chrom <- custom_params$chrom
#   chromstart <- custom_params$chromstart
#   chromend <- custom_params$chromend
#   message("# Reading H3K27ac signals at ", chrom, ":", chromstart, "-", chromend)
#   
#   #
#   bw_dir <- "/Users/chris/Desktop/20240122_H3K27ac_project/output/datahub_UCSC_H3K27ac_0_25m_avg_dvpt/hg38_H3K27ac_0_25m_avg_dvpt"
#   H3K27ac_signal <- list()
#   H3K27ac_timepoint <- paste0(seq(0, 25, 5), "m")
#   for (t in H3K27ac_timepoint) {
#     sample_name <- paste(sep = "_", "H3K27ac", t)
#     message("   > ", sample_name)
#     bw_filename <- paste(sep = "_", "A549", sample_name, "avg_scaled.bw")
#     bw_filepath <- file.path(bw_dir, bw_filename)
#     bw_signal <- readBigwig(file = bw_filepath, params = custom_params)
#     H3K27ac_signal[[sample_name]] <- bw_signal
#   }
#   return(H3K27ac_signal)
# }
# 
# # determine H3K27ac color palette
# H3K27ac_colorpalette <- c("#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B")
# names(H3K27ac_colorpalette) <- paste0("H3K27ac_", seq(0, 25, 5), "m")
# 
# ##### ATAC
# # read ATAC bigwig
# readBigwig_ATAC <- function(custom_params) {
#   chrom <- custom_params$chrom
#   chromstart <- custom_params$chromstart
#   chromend <- custom_params$chromend
#   message("# Reading ATAC signals at : ", chrom, ":", chromstart, "-", chromend)
#   
#   #
#   bw_dir <- "/Users/chris/Desktop/20240122_H3K27ac_project/output/datahub_UCSC_H3K27ac_0_25m_avg_dvpt/hg38_H3K27ac_0_25m_avg_dvpt"
#   ATAC_signal <- list()
#   for (t in c("0h", "1h")) {
#     sample_name <- paste(sep = "_", "ATAC", t)
#     message("   > ", sample_name)
#     bw_filename <- paste(sep = "_", "A549", sample_name, "avg_scaled.bw")
#     bw_filepath <- file.path(bw_dir, bw_filename)
#     bw_signal <- readBigwig(file = bw_filepath, params = custom_params)
#     ATAC_signal[[sample_name]] <- bw_signal
#   }
#   return(ATAC_signal)
# }
# 
# ##### MED1
# # read MED1 bigwig
# readBigwig_MED1 <- function(custom_params) {
#   chrom <- custom_params$chrom
#   chromstart <- custom_params$chromstart
#   chromend <- custom_params$chromend
#   message("# Reading MED1 signals at : ", chrom, ":", chromstart, "-", chromend)
#   
#   #
#   bw_dir <- "/Users/chris/Desktop/20230130_gr_project/output/datahub_UCSC_dvpt/hg38_dvpt"
#   MED1_signal <- list()
#   for (t in c("CTRL", "DEX")) {
#     sample_name <- paste(sep = "_", "MED1", t)
#     message("   > ", sample_name)
#     bw_filename <- paste(sep = "_", "A549", t, "MED1", "pooled_RPKM.bw")
#     bw_filepath <- file.path(bw_dir, bw_filename)
#     bw_signal <- readBigwig(file = bw_filepath, params = custom_params)
#     MED1_signal[[sample_name]] <- bw_signal
#   }
#   return(MED1_signal)
# }
# 
# ##### BRD4
# # read BRD4 bigwig
# readBigwig_BRD4 <- function(custom_params) {
#   chrom <- custom_params$chrom
#   chromstart <- custom_params$chromstart
#   chromend <- custom_params$chromend
#   message("# Reading BRD4 signals at : ", chrom, ":", chromstart, "-", chromend)
#   
#   #
#   bw_dir <- "/Users/chris/Desktop/20230130_gr_project/output/datahub_UCSC_dvpt/hg38_dvpt"
#   BRD4_signal <- list()
#   for (t in c("CTRL", "DEX")) {
#     sample_name <- paste(sep = "_", "BRD4", t)
#     message("   > ", sample_name)
#     bw_filename <- paste(sep = "_", "A549", t, "BRD4", "pooled_RPKM.bw")
#     bw_filepath <- file.path(bw_dir, bw_filename)
#     bw_signal <- readBigwig(file = bw_filepath, params = custom_params)
#     BRD4_signal[[sample_name]] <- bw_signal
#   }
#   return(BRD4_signal)
# }
# 
# ##### GR
# # read GR bigwig
# readBigwig_GR <- function(custom_params) {
#   chrom <- custom_params$chrom
#   chromstart <- custom_params$chromstart
#   chromend <- custom_params$chromend
#   message("# Reading GR signals at : ", chrom, ":", chromstart, "-", chromend)
#   
#   #
#   bw_dir <- "/Users/chris/Desktop/20230130_gr_project/output/datahub_UCSC_dvpt/hg38_dvpt"
#   GR_signal <- list()
#   for (t in c("0h", "1h")) {
#     sample_name <- paste(sep = "_", "GR", t)
#     message("   > ", sample_name)
#     bw_filename <- paste(sep = "_", "A549", "GR", t, "pooled_RPKM.bw")
#     bw_filepath <- file.path(bw_dir, bw_filename)
#     bw_signal <- readBigwig(file = bw_filepath, params = custom_params)
#     GR_signal[[sample_name]] <- bw_signal
#   }
#   return(GR_signal)
# }
# 
# ##### RNA
# # read RNA bigwig
# readBigwig_RNA <- function(custom_params) {
#   chrom <- custom_params$chrom
#   chromstart <- custom_params$chromstart
#   chromend <- custom_params$chromend
#   message("# Reading RNA signals at : ", chrom, ":", chromstart, "-", chromend)
#   
#   #
#   bw_dir <- "/Users/chris/Desktop/20230130_gr_project/output/datahub_UCSC_dvpt/hg38_dvpt"
#   RNA_signal <- list()
#   for (t in c("0h", "0.5h", "1h", "2h", "3h", "4h", "5h", "6h")) {
#     sample_name <- paste(sep = "_", "RNA", t)
#     message("   > ", sample_name)
#     bw_filename <- paste(sep = "_", "A549", "RNA", t, "pooled_RPKM.bw")
#     bw_filepath <- file.path(bw_dir, bw_filename)
#     bw_signal <- readBigwig(file = bw_filepath, params = custom_params)
#     RNA_signal[[sample_name]] <- bw_signal
#   }
#   return(RNA_signal)
# }
# 
# #####
# compute_scalefactor_unit <- function(bw_signal) {
#   # Compute the max score for each df
#   max_scores <- lapply(bw_signal, function(df) max(df$score, na.rm = TRUE))
#   
#   # Identify the max score among the max scores
#   overall_max_score <- max(unlist(max_scores), na.rm = TRUE)
#   
#   # Round the overall max score to the ceiling value
#   rounded_max_score <- ceiling(overall_max_score)
#   
#   # Calculate scale factors for each DataFrame
#   scale_factors <- lapply(max_scores, function(score) rounded_max_score/score)
#   
#   return(scale_factors)
# }
# 
