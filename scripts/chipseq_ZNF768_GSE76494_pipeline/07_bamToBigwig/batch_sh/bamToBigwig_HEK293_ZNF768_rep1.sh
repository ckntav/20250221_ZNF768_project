#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30-ab
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


/cvmfs/soft.mugqic/CentOS6/software/deepTools/deepTools-3.5.4/bin/bamCoverage --extendReads 225 --binSize 10 --smoothLength 30 -p 16 --normalizeUsing RPKM --blackListFileName /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/input/ENCODE_exclusion_list_regions_ENCFF356LFX.bed -b /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE76494-GRCh38_SE/alignment/HEK293_ZNF768_rep1/ZNF768/HEK293_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam -o /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE76494-GRCh38_SE/tracks_byReplicate/HEK293_ZNF768_rep1_RPKM.bw
