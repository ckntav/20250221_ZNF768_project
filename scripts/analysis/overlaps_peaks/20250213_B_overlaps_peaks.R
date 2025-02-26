setwd("/Users/chris/Desktop/20240812_LAPC4_ovERa_project")

library(tidyverse)
library(knitr)
library(kableExtra)
library(viridisLite)
source("scripts/ckn_utils/ckn_utils_load_lapc4_peaks.R")
source("scripts/ckn_utils/ckn_utils_overlaps.R")

#
# Function to compute Jaccard index between two GRanges objects
compute_jaccard_granges <- function(grl) {
  # Ensure input is a GRangesList of length 2
  if (!inherits(grl, "GRangesList") || length(grl) != 2) {
    stop("Input must be a GRangesList with exactly two elements.")
  }
  
  # Extract the two GRanges objects
  gr1 <- grl[[1]]
  gr2 <- grl[[2]]
  
  # Compute intersection (overlapping bases)
  intersection <- sum(width(intersect(gr1, gr2)))
  
  # Compute union (total unique bases covered by either set)
  union <- sum(width(reduce(c(gr1, gr2))))  # Reduce merges overlapping regions
  
  # Compute Jaccard index
  if (union == 0) {
    return(0)  # Avoid division by zero
  } else {
    return(intersection / union)
  }
}

#
antibody_list <- c("AR", "ER")
condition_list <- c("CTRL", "ovERa")
treatment_list <- c("V", "T", "TxD", "D")
reps <- c("rep1", "rep2")

for (antibody in antibody_list) {
  for (condition in condition_list) {
    for (treatment in treatment_list) {
      peaks_i <- load_LAPC4_peaks(antibodies = antibody, conditions = condition, treatments = treatment)
      table_i <- sapply(peaks_i, length) %>% as.data.frame() %>% set_names("number_of_peaks") %>% kbl %>% kable_classic_2(full_width = F, html_font = "Arial", font_size = 20)
      # print(table_i)
      # readline(prompt="Press [enter] to continue")
      
      jaccard_i <- compute_jaccard_granges(peaks_i)
      message("\t> Jaccard index = ", round(jaccard_i, 2))
      
      venn_i <- plotVenn3(peaks_i, labels = TRUE, legend = FALSE, quantities_val = list(type = c("counts")))
      print(venn_i)
      readline(prompt="Press [enter] to continue")
    } 
  }
}
