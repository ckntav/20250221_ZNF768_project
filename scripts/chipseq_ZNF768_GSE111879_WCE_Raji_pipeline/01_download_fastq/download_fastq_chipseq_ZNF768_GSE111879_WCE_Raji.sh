mkdir -p /home/chris11/scratch/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/raw_fastq
cd /home/chris11/scratch/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/raw_fastq

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR684/008/SRR6841438/SRR6841438.fastq.gz -o dl_fastq_SRR6841438.log &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR684/009/SRR6841439/SRR6841439.fastq.gz -o dl_fastq_SRR6841439.log &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR684/000/SRR6841440/SRR6841440.fastq.gz -o dl_fastq_SRR6841440.log &