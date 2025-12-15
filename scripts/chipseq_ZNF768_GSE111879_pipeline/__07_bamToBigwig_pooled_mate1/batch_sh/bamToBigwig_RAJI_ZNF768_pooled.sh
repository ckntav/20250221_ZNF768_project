#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/alignment/RAJI_ZNF768_pooled/ZNF768


samtools merge /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/alignment/RAJI_ZNF768_pooled/ZNF768/RAJI_ZNF768_pooled.ZNF768.sorted.dup.bam /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.bam /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.bam


samtools index /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/alignment/RAJI_ZNF768_pooled/ZNF768/RAJI_ZNF768_pooled.ZNF768.sorted.dup.bam /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/alignment/RAJI_ZNF768_pooled/ZNF768/RAJI_ZNF768_pooled.ZNF768.sorted.dup.bai


bamCoverage --extendReads 225 --binSize 10 --smoothLength 30 -p 16 --samFlagInclude 64 --normalizeUsing RPKM --blackListFileName /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/input/ENCODE_exclusion_list_regions_ENCFF356LFX.bed -b /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/alignment/RAJI_ZNF768_pooled/ZNF768/RAJI_ZNF768_pooled.ZNF768.sorted.dup.bam -o /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/tracks_pooled/RAJI_ZNF768_pooled_RPKM.bw
