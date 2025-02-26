#
today <- get_today()

# retrieve genomic coordinates
genomic_window_bis <- gsub(",", "", genomic_window)
split_window_coord <- str_split(pattern = ":", genomic_window_bis)
chr_i_tmp <- split_window_coord %>% map(1) %>% unlist
# chr_i <- gsub("chr", "", chr_i_tmp)
chr_i <- chr_i_tmp
position_i <- split_window_coord %>% map(2) %>% unlist
start_i <- str_split(pattern = "-", position_i) %>% map(1) %>% unlist %>% as.numeric()
end_i <- str_split(pattern = "-", position_i) %>% map(2) %>% unlist %>% as.numeric()
output_file <- paste(sep = "_", today, i, gene_name)

#
params_i <- pgParams(
  chrom = chr_i, chromstart = start_i, chromend = end_i,
  assembly = hg38_ensembl104,
  # assembly = ensembl104,
  x = 1, just = c("left", "top"),
  width = width_val-1.5, length = width_val-1.5, default.units = "inches"
)
