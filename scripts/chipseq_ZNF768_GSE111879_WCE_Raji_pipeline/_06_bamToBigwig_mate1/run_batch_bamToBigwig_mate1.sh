#!/bin/sh

mkdir -p /home/chris11/projects/def-stbil30/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE/tracks_byReplicate

sbatch scripts/chipseq_ZNF768_GSE111879_WCE_U2OS_pipeline/06_bamToBigwig_mate1/batch_sh/bamToBigwig_mate1_U2OS_WCE_rep1.sh