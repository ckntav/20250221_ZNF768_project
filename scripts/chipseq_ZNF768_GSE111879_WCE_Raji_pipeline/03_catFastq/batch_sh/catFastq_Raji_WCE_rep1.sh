#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


cat /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/raw_fastq/SRR6841438.fastq.gz /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/raw_fastq/SRR6841439.fastq.gz /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/raw_fastq/SRR6841440.fastq.gz > /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/raw_fastq/Raji_WCE_rep1.fastq.gz
