#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/peak_call_custom_0p05/U2OS_ZNF768_rep1/ZNF768


touch /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/peak_call_custom_0p05/U2OS_ZNF768_rep1/ZNF768


macs2 callpeak --format BAMPE --nomodel --gsize 2479938032.8 --treatment /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam -q 0.05 --name /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/peak_call_custom_0p05/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768 >& /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/peak_call_custom_0p05/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.diag.macs.out
