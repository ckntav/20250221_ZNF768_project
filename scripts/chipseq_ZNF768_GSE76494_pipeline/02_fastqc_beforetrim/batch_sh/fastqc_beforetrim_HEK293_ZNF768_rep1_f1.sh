#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30-ab
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE76496-GRCh38_SE/fastqc_beforetrim_output/HEK293_ZNF768_rep1_f1_R1


/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/fastqc/0.12.1/fastqc --outdir /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE76496-GRCh38_SE/fastqc_beforetrim_output/HEK293_ZNF768_rep1_f1_R1 --format fastq /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE76496/raw_fastq/SRR3083355.fastq.gz
