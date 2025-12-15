#!/bin/sh

mkdir -p /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE76494-GRCh38_SE/tracks_byReplicate

sbatch scripts/chipseq_ZNF768_GSE76494_pipeline/07_bamToBigwig/batch_sh/bamToBigwig_HEK293_ZNF768_rep1.sh