setwd("/Users/chris/Desktop/20250221_ZNF768_project")

library(tidyverse)
library(plotgardener)
library(GenomeInfoDb)
source("scripts/ckn_utils/ckn_utils_savePlot.R")
source("scripts/ckn_utils/ckn_utils_ensembl.R")
source("scripts/visualize_genomic_data_4tracks/utils_gviz/gviz_loadRanges.R")

#
output_format <- "tiff"   # "tiff" or "pdf"
tiff_dpi <- 300           # used only when output_format == "tiff"
width_val <- 5
# height_val <- 1.9
height_val <- 3.2
fontsize_val <- 4
height_track <- 0.3

# Load reference genome
ensembl_version <- "v104"
ensembl104 <- get_ensembl(version = "v104")
raw_genes <- genes(ensembl104)
stdchr_genes_tmp <- keepStdChr(raw_genes)
# stdchr_genes <- stdchr_genes_tmp[stdchr_genes_tmp$gene_id %in% active_genes]
stdchr_genes <- stdchr_genes_tmp
seqlevelsStyle(stdchr_genes) <- "UCSC"
# seqlengths_v104 <- seqlengths(stdchr_genes)

# "ZNF768" "ENSG00000169957"
ensembl_gene_list <- c("ENSG00000101412", "ENSG00000007968",
                       "ENSG00000112242", "ENSG00000205250",
                       "ENSG00000133740", "ENSG00000169016",
                       "ENSG00000165891", "ENSG00000129173",
                       "ENSG00000139687", "ENSG00000169957",
                       "ENSG00000099817", "ENSG00000185340", "ENSG00000177455")
symbol_gene_list <- c("E2F1", "E2F2",
                      "E2F3", "E2F4",
                      "E2F5", "E2F6",
                      "E2F7", "E2F8",
                      "RB1", "ZNF768",
                      "POLR2E", "GAS2L1", "CD19")

#
source("scripts/visualize_genomic_data_4tracks/utils_gviz/gviz_visualize_tracks.R")

#
df <- data.frame(ensembl_gene = ensembl_gene_list,
                      symbol = symbol_gene_list) %>% 
  rownames_to_column("rank")


for (i in 1:nrow(df)) {
# for (i in 1:2) {
  ensembl_gene_i <- df %>% dplyr::filter(rank == i) %>% pull(ensembl_gene)
  coord_i <- stdchr_genes[ensembl_gene_i] %>% as.data.frame
  
  kb_around_i <- coord_i$width
  if (kb_around_i > 15000) {
    kb_around_i = 15000
  }

  genomic_window <- paste0(coord_i$seqnames %>% as.character,
                           ":",
                           coord_i$start - kb_around_i,
                           "-",
                           coord_i$end + kb_around_i)
  gene_name <- df %>% dplyr::filter(rank == i) %>% pull(symbol)
  message("# ", i, " | ", gene_name)
  message(" > ", genomic_window)

  #
  source("scripts/visualize_genomic_data_4tracks/utils_gviz/gviz_set_params.R")
  
  #
  source("scripts/visualize_genomic_data_4tracks/utils_gviz/gviz_readBigwig.R")
  
  #
  output_dir <- "output/analysis/tracks_plotgardener_4tracks"
  if (output_format == "tiff") {
    output_filepath <- file.path(output_dir, paste0(output_file, "_", tiff_dpi, "dpi.tiff"))
    tiff(filename = output_filepath, width = width_val, height = height_val,
         units = "in", res = tiff_dpi, compression = "lzw")
  } else if (output_format == "pdf") {
    output_filepath <- file.path(output_dir, paste0(output_file, ".pdf"))
    pdf(file = output_filepath, width = width_val, height = height_val)
  } else {
    stop("Unsupported output_format: ", output_format, " (use 'tiff' or 'pdf')")
  }
  
  # Create pages
  pageCreate(width = width_val, height = height_val, default.units = "inches", showGuides = FALSE)
  
  #
  source("scripts/visualize_genomic_data_4tracks/utils_gviz/gviz_plotSignal.R")
  
  #
  source("scripts/visualize_genomic_data_4tracks/utils_gviz/gviz_plotGenes.R")
  
  #
  source("scripts/visualize_genomic_data_4tracks/utils_gviz/gviz_plotRanges.R")
  
  ## Hide page guides
  # pageGuideHide()
  
  dev.off()
  message(" > Plot (", output_format, ") saved in ", output_filepath)
}
