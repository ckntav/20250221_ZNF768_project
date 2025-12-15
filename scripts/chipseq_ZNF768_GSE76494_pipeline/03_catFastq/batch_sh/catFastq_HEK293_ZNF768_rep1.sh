#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30-ab
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


cat /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE76496/raw_fastq/SRR3083355.fastq.gz /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE76496/raw_fastq/SRR3083356.fastq.gz > /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE76496/raw_fastq/HEK293_ZNF768_rep1.fastq.gz
