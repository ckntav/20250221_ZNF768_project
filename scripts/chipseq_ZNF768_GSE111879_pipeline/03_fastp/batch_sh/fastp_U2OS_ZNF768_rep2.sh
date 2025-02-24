#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/fastp_output/U2OS_ZNF768_rep2


mkdir -p /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/fastp_output


/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/fastp/0.24.0/bin/fastp --in1 /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/raw_fastq/SRR6841442_1.fastq.gz --in2 /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/raw_fastq/SRR6841442_2.fastq.gz --detect_adapter_for_pe --overrepresentation_analysis --overrepresentation_sampling 10 --thread 8 --length_required 48 --out1 /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/fastp_output/U2OS_ZNF768_rep2_1.fastq.gz --out2 /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/fastp_output/U2OS_ZNF768_rep2_2.fastq.gz --html /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/fastp_output/U2OS_ZNF768_rep2/U2OS_ZNF768_rep2_fastp_report.html --json /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/fastp_output/U2OS_ZNF768_rep2/U2OS_ZNF768_rep2_fastp_report.json
