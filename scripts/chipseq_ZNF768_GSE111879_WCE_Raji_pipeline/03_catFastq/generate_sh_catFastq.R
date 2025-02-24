library(tidyverse)

##### module load fastqc

#
fastq_list_filename <- "chipseq_ZNF768_GSE111879_WCE_Raji_fastq_list.txt"
df <- read_tsv(file.path("input", "chipseq_ZNF768_GSE111879_WCE_Raji", fastq_list_filename)) %>% 
  mutate(sample_group = paste(sep = "_", cell_line, target, isogenic_replicate))
fastq_folder <- "chipseq_ZNF768_GSE111879"
output_pipeline_dir <- "chip-pipeline_ZNF768_GSE111879-GRCh38_PE"
script_pipeline_dir <- "chipseq_ZNF768_GSE111879_WCE_Raji_pipeline"
workdir <- "/home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project"

#
header_sh <- c("#!/bin/sh",
               "#SBATCH --time=3:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --ntasks-per-node=1",
               "#SBATCH --cpus-per-task=4",
               "#SBATCH --mem-per-cpu=8G",
               "#SBATCH --account=def-stbil30",
               "#SBATCH --mail-user=christophe.tav@gmail.com",
               "#SBATCH --mail-type=ALL")

sample_name_list <- df$sample_group %>% unique

for (sample_name_i in sample_name_list) {
  # message("##### ", sample_name_i)
  
  df_sample_name <- df %>% dplyr::filter(sample_group == sample_name_i)
  fastq_dir <- file.path(workdir, "raw", fastq_folder, "raw_fastq")
  
  ##### cat R1
  fastq_R1_list <- df_sample_name$fastq_R1_filename
  call_cat_R1_tmp <- c("cat")
  
  for (fastq_R1_filename in fastq_R1_list) {
    fastq_R1_path <- file.path(file.path(workdir, "raw", fastq_folder, "raw_fastq", fastq_R1_filename))
    call_cat_R1_tmp <- paste(call_cat_R1_tmp, fastq_R1_path)
  }
  
  fast_R1_merged_filename <- paste0(sample_name_i, ".fastq.gz")
  fast_R1_merged_filepath <- file.path(file.path(workdir, "raw", fastq_folder, "raw_fastq", fast_R1_merged_filename))
  call_cat_R1 <- paste(call_cat_R1_tmp, ">", fast_R1_merged_filepath)
  
  #
  file_sh1 <- file.path("scripts", script_pipeline_dir , "03_catFastq", "batch_sh",
                        paste0("catFastq_", sample_name_i, ".sh"))
  message("sbatch ", file_sh1)
  fileConn1 <- file(file_sh1)
  writeLines(c(header_sh, "\n", call_cat_R1), fileConn1)
  close(fileConn1)
  
  ##### cat R2
  # fastq_R2_list <- df_sample_name$fastq_R2_filename
  # call_cat_R2_tmp <- c("cat")
  # 
  # for (fastq_R2_filename in fastq_R2_list) {
  #   fastq_R2_path <- file.path(file.path(workdir, "raw", fastq_folder, "raw_fastq", fastq_R2_filename))
  #   call_cat_R2_tmp <- paste(call_cat_R2_tmp, fastq_R2_path)
  # }
  # 
  # fast_R2_merged_filename <- paste0(sample_name_i, "_R2.fastq.gz")
  # fast_R2_merged_filepath <- file.path(file.path(workdir, "raw", fastq_folder, "raw_fastq", fast_R2_merged_filename))
  # call_cat_R2 <- paste(call_cat_R2_tmp, ">", fast_R2_merged_filepath)
  # 
  # #
  # file_sh2 <- file.path("scripts", script_pipeline_dir , "03_catFastq", "batch_sh",
  #                       paste0("catFastq_", sample_name_i, "_R2.sh"))
  # message("sbatch ", file_sh2)
  # fileConn2 <- file(file_sh2)
  # writeLines(c(header_sh, "\n", call_cat_R2), fileConn2)
  # close(fileConn2)
}
