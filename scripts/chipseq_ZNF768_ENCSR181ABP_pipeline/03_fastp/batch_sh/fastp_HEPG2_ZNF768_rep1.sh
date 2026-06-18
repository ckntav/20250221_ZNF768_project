#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-adroit
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_ENCSR181ABP-GRCh38_SE/fastp_output/HEPG2_ZNF768_rep1


mkdir -p /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_ENCSR181ABP/fastp_output


/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/fastp/0.24.0/bin/fastp --in1 /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_ENCSR181ABP/raw_fastq/ENCFF300JRW.fastq.gz --overrepresentation_analysis --overrepresentation_sampling 10 --thread 8 --length_required 74 --out1 /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_ENCSR181ABP/fastp_output/HEPG2_ZNF768_rep1_1.fastq.gz --html /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_ENCSR181ABP-GRCh38_SE/fastp_output/HEPG2_ZNF768_rep1/HEPG2_ZNF768_rep1_fastp_report.html --json /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_ENCSR181ABP-GRCh38_SE/fastp_output/HEPG2_ZNF768_rep1/HEPG2_ZNF768_rep1_fastp_report.json
