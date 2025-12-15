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
# all_HEK293_peaks <-all_peaks[grepl("U2OS", names(all_peaks))] # to_change