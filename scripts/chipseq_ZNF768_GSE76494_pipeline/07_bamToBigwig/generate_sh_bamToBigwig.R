setwd("/Users/chris/Desktop/20250221_ZNF768_project")

library(tidyverse)

##### module load fastp/0.23.4

#
pi_id <- "def-stbil30-ab"

#
fastq_list_filename <- "chipseq_ZNF768_GSE76496_fastq_list.txt"
df <- read_tsv(file.path("input", "chipseq_ZNF768_GSE76496", fastq_list_filename))
output_pipeline_dir <- "chip-pipeline_ZNF768_GSE76494-GRCh38_SE"
script_pipeline_dir <- "chipseq_ZNF768_GSE76494_pipeline"
workdir <- "/home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project"

alignment_dir <- file.path("output", output_pipeline_dir, "alignment")
tracks_dir <- file.path("output", output_pipeline_dir, "tracks_byReplicate")
message("mkdir -p ", file.path(workdir, tracks_dir))


header_sh <- c("#!/bin/sh",
               "#SBATCH --time=3:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --ntasks-per-node=1",
               "#SBATCH --cpus-per-task=16",
               "#SBATCH --mem-per-cpu=8G",
               paste0("#SBATCH --account=", pi_id),
               "#SBATCH --mail-user=christophe.tav@gmail.com",
               "#SBATCH --mail-type=ALL")

ENCODE_blacklist <- "input/ENCODE_exclusion_list_regions_ENCFF356LFX.bed"
bamCoverage_path <- "/cvmfs/soft.mugqic/CentOS6/software/deepTools/deepTools-3.5.4/bin/bamCoverage"
norm <- "RPKM"

if (TRUE) {
  sample_name <- "HEK293_ZNF768_rep1"
  antibody_factor <- "ZNF768"
  # message(" > ", sample_name)
  
  bam <- file.path(workdir, alignment_dir, sample_name, antibody_factor, paste0(sample_name, ".",  antibody_factor, ".sorted.dup.filtered.bam"))
  bw <- file.path(workdir, tracks_dir, paste0(sample_name, "_", norm, ".bw"))
  # message(paste0(sample_name, "_", norm, "_1.bw"))
  
  # print(paste0(" > bam : ", bam))
  # print(paste0(" > bigwig : ", bw))
  # message(basename(bw))
  call_bamCoverage <- paste(bamCoverage_path, "--extendReads", 225, "--binSize", 10,
                            "--smoothLength", 30,
                            "-p", 16,
                            # "--samFlagInclude", 64,
                            "--normalizeUsing", norm,
                            "--blackListFileName", file.path(workdir, ENCODE_blacklist),
                            "-b", bam,
                            "-o", bw)
  # message(call_bamCoverage)
  # system(call_bamCoverage)
  
  file_sh <- file.path("scripts", script_pipeline_dir, "07_bamToBigwig", "batch_sh",
                       paste0("bamToBigwig_", sample_name, ".sh"))
  message("sbatch ", file_sh)
  fileConn <- file(file_sh)
  writeLines(c(header_sh, "\n", call_bamCoverage), fileConn)
  close(fileConn)
}