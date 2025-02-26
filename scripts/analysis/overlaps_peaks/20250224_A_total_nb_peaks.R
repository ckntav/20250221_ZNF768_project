library(tidyverse)
library(knitr)
library(kableExtra)
source("scripts/ckn_utils/ckn_utils_load_peaks.R")
source("scripts/ckn_utils/ckn_utils_overlaps.R")

#
all_peaks <- load_ZNF768_peaks()

sapply(all_peaks, length) %>%
  as.data.frame %>%
  set_names("number_of_peaks") %>%
  kbl %>%
  kable_paper(full_width = F)
