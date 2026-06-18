source("scripts/ckn_utils/ckn_utils_load_peaks.R")

all_peaks <- load_ZNF768_peaks()

# all_H3_macs2peaks <- load_h3_macs2peaks()
# all_H3_gopeaks <- load_h3_gopeaks()
# H3K4me3_macs2peaks <- all_H3_macs2peaks[grepl("H3K4me3", names(all_H3_macs2peaks))]
# names(H3K4me3_macs2peaks) <- paste0(names(H3K4me3_macs2peaks), "_macs2")
# H3K27me3_macs2peaks <- all_H3_macs2peaks[grepl("H3K27me3", names(all_H3_macs2peaks))]
# names(H3K27me3_macs2peaks) <- paste0(names(H3K27me3_macs2peaks), "_macs2")
# H3K4me3_gopeaks <- all_H3_gopeaks[grepl("H3K4me3", names(all_H3_gopeaks))]
# H3K27me3_gopeaks <- all_H3_gopeaks[grepl("H3K27me3", names(all_H3_gopeaks))]

all_RAJI_peaks <- all_peaks[grepl("RAJI", names(all_peaks))]
all_U2OS_peaks <-all_peaks[grepl("U2OS", names(all_peaks))]

#
all_HEK293_peaks <- GRangesList()
HEK293_ZNF768_rep1_peaks <- rtracklayer::import("output/GEO_peak_GSE76496/HEK293_ZNF768_rep1.sorted.bed")
all_HEK293_peaks[["HEK293_ZNF768_rep1"]] <- HEK293_ZNF768_rep1_peaks

#
read_narrowpeak <- function(path) {
  df <- data.table::fread(path)
  colnames(df) <- c("chr", "start", "end", "name", "score", "strand",
                    "signalValue", "pValue", "qValue", "peak")
  GRanges(
    seqnames    = df$chr,
    ranges      = IRanges(start = df$start + 1L, end = df$end),
    strand      = "*",           # strand column is "." so force *
    name        = df$name,
    score       = as.integer(df$score),
    signalValue = df$signalValue,
    pValue      = df$pValue,
    qValue      = df$qValue,
    peak        = df$peak
  )
}

all_HEPG2_peaks <- GRangesList()
HEPG2_ZNF768_rep1_peaks <- read_narrowpeak("output/ENCODE_peak_ENCSR181ABP/HEPG2_ZNF768_rep1_peaks_ENCFF597ISS.bed.gz")
HEPG2_ZNF768_rep2_peaks <- read_narrowpeak("output/ENCODE_peak_ENCSR181ABP/HEPG2_ZNF768_rep2_peaks_ENCFF675UFF.bed.gz")
all_HEPG2_peaks[["HEPG2_ZNF768_rep1"]] <- HEPG2_ZNF768_rep1_peaks
all_HEPG2_peaks[["HEPG2_ZNF768_rep2"]] <- HEPG2_ZNF768_rep2_peaks
