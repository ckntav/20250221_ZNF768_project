#!/bin/sh

mkdir -p /home/chris11/projects/def-stbil30-ab/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_ENCSR181ABP-GRCh38_SE/tracks_pooled

sbatch scripts/chipseq_ZNF768_ENCSR181ABP_pipeline/07_bamToBigwig_pooled/batch_sh/bamToBigwig_HEPG2_ZNF768_pooled.sh