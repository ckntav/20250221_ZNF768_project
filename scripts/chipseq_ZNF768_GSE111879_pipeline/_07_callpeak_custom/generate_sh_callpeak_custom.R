library(tidyverse)

##### mugqic
##### module load mugqic/python/3.7.3 mugqic/MACS2/2.2.7.1

#
fastq_list_filename <- "chipseq_ZNF768_GSE111879_fastq_list.txt"
df <- read_tsv(file.path("input", "chipseq_ZNF768_GSE111879", fastq_list_filename)) %>% 
  mutate(antibody_factor = ifelse(target == "WCE", "Input", target)) %>% 
  dplyr::filter(target != "WCE")
output_pipeline_dir <- "chip-pipeline_ZNF768_GSE111879-GRCh38_PE"
script_pipeline_dir <- "chipseq_ZNF768_GSE111879_pipeline"
workdir <- "/home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project"

# alignment_dir <- "output/chip-pipeline-GRCh38_1rep/alignment"
output_dir <- file.path("output", output_pipeline_dir)

header_sh <- c("#!/bin/sh",
               "#SBATCH --time=3:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --ntasks-per-node=1",
               "#SBATCH --cpus-per-task=16",
               "#SBATCH --mem-per-cpu=8G",
               "#SBATCH --account=def-stbil30",
               "#SBATCH --mail-user=christophe.tav@gmail.com",
               "#SBATCH --mail-type=ALL")

qval <- 0.05
qval_chr <- gsub("\\.", "p", qval)
format_aln <- "BAMPE"

for (i in 1:nrow(df)) {
  sample_name <- df$sample_name[i]
  cell_line_i <- df$cell_line[i]
  antibody_factor <- df$antibody_factor[i]
  # condition_i <- df$condition[i]
  # treatment_i <- df$treatment[i]
  rep_i <- df$isogenic_replicate[i]
  
  output_alnrep_dir <- file.path(output_dir, "alignment")
  
  basename <- paste(cell_line_i, antibody_factor, rep_i, sep = "_")
  # message(" > ", basename)
  
  peak_dir <- file.path(workdir, output_dir, paste(sep = "_", "peak_call_custom", qval_chr), basename, antibody_factor)
  
  bam_trt <- file.path(workdir, output_alnrep_dir, basename, antibody_factor,
                       paste0(basename, ".",  antibody_factor, ".sorted.dup.filtered.bam"))
  
  call_mkdir <- paste("mkdir", "-p", peak_dir)
  call_touch <- paste("touch", peak_dir)
  
  call_callpeak <- paste("macs2", "callpeak",
                         "--format", format_aln,
                         "--nomodel",
                         "--gsize 2479938032.8",
                         "--treatment", bam_trt,
                         # "--control", bam_ctrl,
                         # "--nolambda",
                         "-q", qval,
                         "--name", paste0(peak_dir, "/", basename, ".", antibody_factor),
                         ">&", paste0(peak_dir, "/", basename, ".", antibody_factor, ".diag.macs.out"))
  
  # message(call_callpeak)
  
  file_sh <- file.path("scripts", script_pipeline_dir, "07_callpeak_custom", "batch_sh",
                       paste0("callpeak_custom_", qval_chr, "_", basename, ".sh"))
  message("sbatch ", file_sh)
  fileConn <- file(file_sh)
  writeLines(c(header_sh, "\n", call_mkdir, "\n", call_touch, "\n", call_callpeak), fileConn)
  close(fileConn)
}
