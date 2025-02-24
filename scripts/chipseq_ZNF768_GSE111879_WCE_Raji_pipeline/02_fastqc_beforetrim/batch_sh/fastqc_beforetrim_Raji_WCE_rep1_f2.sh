#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/fastqc_beforetrim_output/Raji_WCE_rep1_f2_R1


/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/fastqc/0.12.1/fastqc --outdir /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/fastqc_beforetrim_output/Raji_WCE_rep1_f2_R1 --format fastq /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111880/raw_fastq/SRR6841439.fastq.gz
