mkdir -p $SCRATCH/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE

genpipes chipseq --job-scheduler slurm -s 1-13 \
    --log debug \
    --readsets raw/chipseq_ZNF768_GSE111879/readset_chipseq_ZNF768_GSE111879_20250221.txt \
    --output-dir output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE \
    --no-json \
    --config $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.base.ini \
        $MUGQIC_PIPELINES_HOME/pipelines/common_ini/narval.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini \
    --genpipes_file scripts/chipseq_ZNF768_GSE111879_pipeline/05_mugqic/20250221_chipseq_ZNF768_GSE111879_GRCh38_pipeline.sh