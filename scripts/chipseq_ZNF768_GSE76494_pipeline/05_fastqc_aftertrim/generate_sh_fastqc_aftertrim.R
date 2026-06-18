library(tidyverse)

##### module load fastqc

#
pi_id <- "def-adroit"

#
fastq_list_filename <- "chipseq_ZNF768_ENCSR181ABP_fastq_list.txt"
df <- read_tsv(file.path("input", "chipseq_ZNF768_ENCSR181ABP", fastq_list_filename))
output_pipeline_dir <- "chip-pipeline_ZNF768_ENCSR181ABP-GRCh38_SE"
script_pipeline_dir <- "chipseq_ZNF768_ENCSR181ABP_pipeline"
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

