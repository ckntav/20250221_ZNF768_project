mkdir -p $SCRATCH/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE76494-GRCh38_SE

genpipes chipseq --job-scheduler slurm -s 1-13 \
    --log debug \
    --readsets raw/chipseq_ZNF768_GSE76496/readset_chipseq_ZNF768_GSE76494_20251213.txt \
    --output-dir output/chip-pipeline_ZNF768_GSE76494-GRCh38_SE \
    --config $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.base.ini \
        $MUGQIC_PIPELINES_HOME/pipelines/common_ini/narval.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini \
    --genpipes_file scripts/chipseq_ZNF768_GSE76494_pipeline/06_mugqic/20251213_chipseq_ZNF768_GSE76494_GRCh38_pipeline.sh