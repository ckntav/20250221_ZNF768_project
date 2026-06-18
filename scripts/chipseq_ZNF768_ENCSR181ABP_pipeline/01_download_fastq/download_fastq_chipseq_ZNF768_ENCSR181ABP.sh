mkdir -p /home/chris11/scratch/20250221_ZNF768_project/raw/chipseq_ZNF768_ENCSR181ABP/raw_fastq
cd /home/chris11/scratch/20250221_ZNF768_project/raw/chipseq_ZNF768_ENCSR181ABP/raw_fastq

wget https://www.encodeproject.org/files/ENCFF300JRW/@@download/ENCFF300JRW.fastq.gz -o dl_fastq_ENCFF300JRW.log &
wget https://www.encodeproject.org/files/ENCFF020PJB/@@download/ENCFF020PJB.fastq.gz -o dl_fastq_ENCFF020PJB.log &