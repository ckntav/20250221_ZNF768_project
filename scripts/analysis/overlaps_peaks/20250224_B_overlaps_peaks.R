library(tidyverse)
library(knitr)
library(kableExtra)
library(viridisLite)
source("scripts/ckn_utils/ckn_utils_overlaps.R")
source("scripts/ckn_utils/ckn_utils_load_peaks.R")

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
all_peaks <- load_ZNF768_peaks()
viridis4 <- viridis(n = 4)
plotVenn3(all_peaks, labels = FALSE, legend = TRUE, quantities_val = list(type = c("counts")))
combMat_all_peaks <- generate_comb_mat(all_peaks)
displayUpSet(combMat_all_peaks)

# RAJI
raji_peaks <- all_peaks[grepl("RAJI", names(all_peaks))]
names(raji_peaks)
sapply(raji_peaks, length) %>%
  as.data.frame %>%
  set_names("number_of_peaks") %>%
  kbl %>%
  kable_paper(full_width = F)
jaccard_raji <- compute_jaccard_granges(raji_peaks)
message("\t> Jaccard index = ", round(jaccard_raji, 2))
venn_raji <- plotVenn3(raji_peaks, labels = TRUE, legend = FALSE, quantities_val = list(type = c("counts")))
print(venn_raji)

# U2OS
u2os_peaks <- all_peaks[grepl("U2OS", names(all_peaks))]
names(u2os_peaks)
sapply(u2os_peaks, length) %>%
  as.data.frame %>%
  set_names("number_of_peaks") %>%
  kbl %>%
  kable_paper(full_width = F)
jaccard_u2os <- compute_jaccard_granges(u2os_peaks)
message("\t> Jaccard index = ", round(jaccard_u2os, 2))
venn_u2os <- plotVenn3(u2os_peaks, labels = TRUE, legend = FALSE, quantities_val = list(type = c("counts")))
print(venn_u2os)
