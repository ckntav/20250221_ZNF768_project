library(tidyverse)

##### module load fastqc

#
pi_id <- "def-stbil30-ab"

#
fastq_list_filename <- "chipseq_ZNF768_GSE76496_fastq_list.txt"
df <- read_tsv(file.path("input", "chipseq_ZNF768_GSE76496", fastq_list_filename))
output_pipeline_dir <- "chip-pipeline_ZNF768_GSE76496-GRCh38_SE"
script_pipeline_dir <- "chipseq_ZNF768_GSE76494_pipeline"
workdir <- "/home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project"

#
header_sh <- c("#!/bin/sh",
               "#SBATCH --time=3:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --ntasks-per-node=1",
               "#SBATCH --cpus-per-task=4",
               "#SBATCH --mem-per-cpu=8G",
               paste0("#SBATCH --account=", pi_id),
               "#SBATCH --mail-user=christophe.tav@gmail.com",
               "#SBATCH --mail-type=ALL")

fastqc_path <- "/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/fastqc/0.12.1/fastqc"

if (TRUE) {
  sample_name <- "HEK293_ZNF768_rep1"
  fastq_folder_i <- df$fastq_folder %>% unique
  
  output_fastp_dir <- file.path(workdir, "raw", fastq_folder_i, "fastp_output")
  
  ##### R1
  # message(" > R1")
  fastq_R1_filename <- paste0(sample_name, ".fastq.gz")
  fastq_R1_filepath <- file.path(output_fastp_dir, fastq_R1_filename )
  output_fastqc_R1_path <- file.path(workdir, "output", output_pipeline_dir, "fastqc_aftertrim_output", paste(sep = "_", sample_name, "R1"))
  call_mkdir_R1 <- paste("mkdir", "-p", output_fastqc_R1_path)
  
  #
  call_fastqc_R1 <- paste(fastqc_path,
                          "--outdir", output_fastqc_R1_path,
                          "--format", "fastq",
                          fastq_R1_filepath)
  
  #
  file_sh1 <- file.path("scripts", script_pipeline_dir, "05_fastqc_aftertrim/batch_sh",
                        paste0("fastqc_aftertrim_", sample_name, "_R1.sh"))
  message("sbatch ", file_sh1)
  fileConn1 <- file(file_sh1)
  writeLines(c(header_sh, "\n", call_mkdir_R1, "\n", call_fastqc_R1), fileConn1)
  close(fileConn1)
}

