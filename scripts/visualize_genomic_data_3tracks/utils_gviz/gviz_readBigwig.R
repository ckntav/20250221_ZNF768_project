#
RAJI_bw <- readBigwig_RAJI(custom_params = params_i)
RAJI_maxRange <- set_maxRange(RAJI_bw)
if (RAJI_maxRange == 0) {
  RAJI_maxRange = 1
}

#
U2OS_bw <- readBigwig_U2OS(custom_params = params_i)
U2OS_maxRange <- set_maxRange(U2OS_bw)
if (U2OS_maxRange == 0) {
  U2OS_maxRange = 1
}

#
HEK293_bw <- readBigwig_HEK293(custom_params = params_i)
HEK293_maxRange <- set_maxRange(HEK293_bw)
if (HEK293_maxRange == 0) {
  HEK293_maxRange = 1
}

#
# RNA_bw <- readBigwig_RNA(custom_params = params_i)
# RNA_scalefactor <- compute_scalefactor(RNA_bw)

# #
# ATAC_bw <- readBigwig_ATAC(custom_params = params_i)
# ATAC_scalefactor <- compute_scalefactor_unit(ATAC_bw)
# #
# MED1_bw <- readBigwig_MED1(custom_params = params_i)
# MED1_scalefactor <- compute_scalefactor(MED1_bw)
# #
# BRD4_bw <- readBigwig_BRD4(custom_params = params_i)
# BRD4_scalefactor <- compute_scalefactor(BRD4_bw)
# #
# GR_bw <- readBigwig_GR(custom_params = params_i)
# GR_scalefactor <- compute_scalefactor(GR_bw)