#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-adroit
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_ENCSR181ABP-GRCh38_SE/alignment/HEPG2_ZNF768_pooled/ZNF768


samtools merge /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_ENCSR181ABP-GRCh38_SE/alignment/HEPG2_ZNF768_pooled/ZNF768/HEPG2_ZNF768_pooled.ZNF768.sorted.dup.bam /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_ENCSR181ABP-GRCh38_SE/alignment/HEPG2_ZNF768_rep1/ZNF768/HEPG2_ZNF768_rep1.ZNF768.sorted.dup.bam /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_ENCSR181ABP-GRCh38_SE/alignment/HEPG2_ZNF768_rep2/ZNF768/HEPG2_ZNF768_rep2.ZNF768.sorted.dup.bam


samtools index /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_ENCSR181ABP-GRCh38_SE/alignment/HEPG2_ZNF768_pooled/ZNF768/HEPG2_ZNF768_pooled.ZNF768.sorted.dup.bam /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_ENCSR181ABP-GRCh38_SE/alignment/HEPG2_ZNF768_pooled/ZNF768/HEPG2_ZNF768_pooled.ZNF768.sorted.dup.bai


bamCoverage --extendReads 225 --binSize 10 --smoothLength 30 -p 16 --normalizeUsing RPKM --blackListFileName /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/input/ENCODE_exclusion_list_regions_ENCFF356LFX.bed -b /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_ENCSR181ABP-GRCh38_SE/alignment/HEPG2_ZNF768_pooled/ZNF768/HEPG2_ZNF768_pooled.ZNF768.sorted.dup.bam -o /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_ENCSR181ABP-GRCh38_SE/tracks_pooled/HEPG2_ZNF768_pooled_RPKM.bw
