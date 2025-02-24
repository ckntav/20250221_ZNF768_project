library(tidyverse)

##### module load fastp/0.23.4

#
fastq_list_filename <- "chipseq_ZNF768_GSE111879_fastq_list.txt"
df <- read_tsv(file.path("input", "chipseq_ZNF768_GSE111879", fastq_list_filename)) %>% 
  mutate(antibody_factor = ifelse(target == "WCE", "Input", target)) %>% 
  mutate(sample_group = paste(sep = "_", cell_line, target))
output_pipeline_dir <- "chip-pipeline_ZNF768_GSE111879-GRCh38_PE"
script_pipeline_dir <- "chipseq_ZNF768_GSE111879_pipeline"
workdir <- "/home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project"

alignment_dir <- file.path("output", output_pipeline_dir, "alignment")
tracks_dir <- file.path("output", output_pipeline_dir, "tracks_pooled")
message("mkdir -p ", file.path(workdir, tracks_dir))


header_sh <- c("#!/bin/sh",
               "#SBATCH --time=3:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --ntasks-per-node=1",
               "#SBATCH --cpus-per-task=16",
               "#SBATCH --mem-per-cpu=8G",
               "#SBATCH --account=def-stbil30",
               "#SBATCH --mail-user=christophe.tav@gmail.com",
               "#SBATCH --mail-type=ALL")

ENCODE_blacklist <- "input/ENCODE_exclusion_list_regions_ENCFF356LFX.bed"

norm <- "RPKM"

#
sample_group_list <- df %>% pull(sample_group) %>% unique

for (sample_group in sample_group_list) {
  # message("# ", sample_group)
  basename_merge <- paste(sample_group, "pooled", sep = "_")
  
  outBam_folder_filepath <- file.path(workdir, alignment_dir, basename_merge, "ZNF768")
  call_mkdir <- paste("mkdir", "-p", outBam_folder_filepath)
  # message(call_mkdir)
  # message()
  
  outBam_merge <- file.path(outBam_folder_filepath, paste0(basename_merge, ".", "ZNF768", ".sorted.dup.bam"))
  
  # message("  ### Pooled bam between replicates")
  # message("    > outBam : ", outBam_merge)
  call_samtools_merge <- paste("samtools", "merge", outBam_merge)
  
  for (rep in c("rep1", "rep2")) {
    basename <- paste(sample_group, rep, sep = "_")
    bam <- file.path(workdir, alignment_dir, basename, "ZNF768", paste0(basename, ".", "ZNF768", ".sorted.dup.bam"))
    # message("    > bam ", rep, " : ", bam)
    
    call_samtools_merge <- paste(call_samtools_merge, bam)
  }
  
  # message(call_samtools_merge)
  # message()
  
  ##### Index pooled bam
  # message("  ### Index pooled bam")
  outBam_index <- gsub(".bam", ".bai", outBam_merge)
  # message("    > outBam pooled : ", outBam_merge)
  # message("    > outBam index : ", outBam_index)
  call_samtools_index <- paste("samtools", "index", outBam_merge, outBam_index)
  # message(call_samtools_index)
  # message()
  
  
  ##### bamToBigWig
  # message("  ### bamToBigWig")
  bam <- outBam_merge
  bw <- file.path(workdir, tracks_dir, paste0(basename_merge, "_", norm, ".bw"))
  # message("    > bam : ", bam)
  # message("    > bigwig : ", bw)
  # 
  call_bamCoverage <- paste("bamCoverage", "--extendReads", 225, "--binSize", 10,
                            "--smoothLength", 30,
                            "-p", 16,
                            "--samFlagInclude", 64,
                            "--normalizeUsing", norm,
                            "--blackListFileName", file.path(workdir, ENCODE_blacklist),
                            "-b", bam,
                            "-o", bw)
  # message(call_bamCoverage)
  # message()
  # system(call_bamCoverage)
  
  ##### Generate script
  file_sh <- file.path("scripts", script_pipeline_dir, "07_bamToBigwig_pooled_mate1", "batch_sh",
                       paste0("bamToBigwig_", basename_merge, ".sh"))
  message("sbatch ", file_sh)
  fileConn <- file(file_sh)
  writeLines(c(header_sh, "\n",
               call_mkdir, "\n",
               call_samtools_merge, "\n",
               call_samtools_index, "\n",
               call_bamCoverage),
             fileConn)
  close(fileConn)
}