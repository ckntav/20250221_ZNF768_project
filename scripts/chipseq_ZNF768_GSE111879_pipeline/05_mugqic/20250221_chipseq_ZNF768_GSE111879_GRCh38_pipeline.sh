#!/bin/bash
# Exit immediately on error

set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq SLURM Job Submission Bash script
# Version: 5.1.0
# Created on: 2025-02-21T16.40.33
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 4 jobs
#   merge_trimmomatic_stats: 1 job
#   mapping_bwa_mem_sambamba: 4 jobs
#   sambamba_merge_bam_files: 4 jobs
#   sambamba_mark_duplicates: 4 jobs
#   sambamba_view_filter: 4 jobs
#   bedtools_blacklist_filter: 0 job... skipping
#   metrics: 9 jobs
#   homer_make_tag_directory: 4 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 9 jobs
#   macs2_callpeak: 8 jobs
#   TOTAL: 52 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/lustre06/project/6001942/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=2025-02-21T16.40.33
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq.chipseq.job_list.$TIMESTAMP
export CONFIG_FILES="/cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-5.1.0/lib/python3.12/site-packages/genpipes/pipelines/chipseq/chipseq.base.ini,/cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-5.1.0/lib/python3.12/site-packages/genpipes/pipelines/common_ini/narval.ini,/cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR
#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.ChIP_seq_RAJI_ZNF768_rep1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ChIP_seq_RAJI_ZNF768_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ChIP_seq_RAJI_ZNF768_rep1.c8e76ca1cf8fd9a0567c38aa256c21d2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.ChIP_seq_RAJI_ZNF768_rep1.c8e76ca1cf8fd9a0567c38aa256c21d2.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.39 && \
mkdir -p trim/RAJI_ZNF768_rep1/ZNF768 && \
touch trim/RAJI_ZNF768_rep1/ZNF768 && \
mkdir -p metrics/multiqc_inputs && \
touch metrics/multiqc_inputs && \
`cat > trim/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1.trim.adapters.fa << END
>Prefix/1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>Prefix/2
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
END
` && \
java -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx20G -jar $TRIMMOMATIC_JAR PE \
  -threads 8 \
  -phred33 \
  /lustre06/project/6001942/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/fastp_output/RAJI_ZNF768_rep1_1.fastq.gz \
  /lustre06/project/6001942/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/fastp_output/RAJI_ZNF768_rep1_2.fastq.gz \
  trim/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1.trim.pair1.fastq.gz \
  trim/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1.trim.single1.fastq.gz \
  trim/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1.trim.pair2.fastq.gz \
  trim/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1.trim.single2.fastq.gz \
  ILLUMINACLIP:trim/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1.trim.log && \
ln -s -f \
  ../../trim/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1.trim.log \
  metrics/multiqc_inputs/ChIP_seq_RAJI_ZNF768_rep1.trim.log
trimmomatic.ChIP_seq_RAJI_ZNF768_rep1.c8e76ca1cf8fd9a0567c38aa256c21d2.mugqic.done
chmod 755 $COMMAND
trimmomatic_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 20G -c 8 -N 1    | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trimmomatic_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.ChIP_seq_RAJI_ZNF768_rep2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ChIP_seq_RAJI_ZNF768_rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ChIP_seq_RAJI_ZNF768_rep2.444663bcbf4f6be773373883833572e2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.ChIP_seq_RAJI_ZNF768_rep2.444663bcbf4f6be773373883833572e2.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.39 && \
mkdir -p trim/RAJI_ZNF768_rep2/ZNF768 && \
touch trim/RAJI_ZNF768_rep2/ZNF768 && \
mkdir -p metrics/multiqc_inputs && \
touch metrics/multiqc_inputs && \
`cat > trim/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2.trim.adapters.fa << END
>Prefix/1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>Prefix/2
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
END
` && \
java -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx20G -jar $TRIMMOMATIC_JAR PE \
  -threads 8 \
  -phred33 \
  /lustre06/project/6001942/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/fastp_output/RAJI_ZNF768_rep2_1.fastq.gz \
  /lustre06/project/6001942/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/fastp_output/RAJI_ZNF768_rep2_2.fastq.gz \
  trim/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2.trim.pair1.fastq.gz \
  trim/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2.trim.single1.fastq.gz \
  trim/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2.trim.pair2.fastq.gz \
  trim/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2.trim.single2.fastq.gz \
  ILLUMINACLIP:trim/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2.trim.log && \
ln -s -f \
  ../../trim/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2.trim.log \
  metrics/multiqc_inputs/ChIP_seq_RAJI_ZNF768_rep2.trim.log
trimmomatic.ChIP_seq_RAJI_ZNF768_rep2.444663bcbf4f6be773373883833572e2.mugqic.done
chmod 755 $COMMAND
trimmomatic_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 20G -c 8 -N 1    | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trimmomatic_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.ChIP_seq_U2OS_ZNF768_rep1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ChIP_seq_U2OS_ZNF768_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ChIP_seq_U2OS_ZNF768_rep1.66721c8efa305eecea5c78eea372b30f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.ChIP_seq_U2OS_ZNF768_rep1.66721c8efa305eecea5c78eea372b30f.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.39 && \
mkdir -p trim/U2OS_ZNF768_rep1/ZNF768 && \
touch trim/U2OS_ZNF768_rep1/ZNF768 && \
mkdir -p metrics/multiqc_inputs && \
touch metrics/multiqc_inputs && \
`cat > trim/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1.trim.adapters.fa << END
>Prefix/1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>Prefix/2
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
END
` && \
java -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx20G -jar $TRIMMOMATIC_JAR PE \
  -threads 8 \
  -phred33 \
  /lustre06/project/6001942/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/fastp_output/U2OS_ZNF768_rep1_1.fastq.gz \
  /lustre06/project/6001942/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/fastp_output/U2OS_ZNF768_rep1_2.fastq.gz \
  trim/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1.trim.pair1.fastq.gz \
  trim/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1.trim.single1.fastq.gz \
  trim/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1.trim.pair2.fastq.gz \
  trim/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1.trim.single2.fastq.gz \
  ILLUMINACLIP:trim/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1.trim.log && \
ln -s -f \
  ../../trim/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1.trim.log \
  metrics/multiqc_inputs/ChIP_seq_U2OS_ZNF768_rep1.trim.log
trimmomatic.ChIP_seq_U2OS_ZNF768_rep1.66721c8efa305eecea5c78eea372b30f.mugqic.done
chmod 755 $COMMAND
trimmomatic_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 20G -c 8 -N 1    | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trimmomatic_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.ChIP_seq_U2OS_ZNF768_rep2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ChIP_seq_U2OS_ZNF768_rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ChIP_seq_U2OS_ZNF768_rep2.eca5d0f1e9d90118c2dbb58a91b94414.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.ChIP_seq_U2OS_ZNF768_rep2.eca5d0f1e9d90118c2dbb58a91b94414.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.39 && \
mkdir -p trim/U2OS_ZNF768_rep2/ZNF768 && \
touch trim/U2OS_ZNF768_rep2/ZNF768 && \
mkdir -p metrics/multiqc_inputs && \
touch metrics/multiqc_inputs && \
`cat > trim/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2.trim.adapters.fa << END
>Prefix/1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>Prefix/2
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
END
` && \
java -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx20G -jar $TRIMMOMATIC_JAR PE \
  -threads 8 \
  -phred33 \
  /lustre06/project/6001942/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/fastp_output/U2OS_ZNF768_rep2_1.fastq.gz \
  /lustre06/project/6001942/chris11/20250221_ZNF768_project/raw/chipseq_ZNF768_GSE111879/fastp_output/U2OS_ZNF768_rep2_2.fastq.gz \
  trim/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2.trim.pair1.fastq.gz \
  trim/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2.trim.single1.fastq.gz \
  trim/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2.trim.pair2.fastq.gz \
  trim/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2.trim.single2.fastq.gz \
  ILLUMINACLIP:trim/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2.trim.log && \
ln -s -f \
  ../../trim/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2.trim.log \
  metrics/multiqc_inputs/ChIP_seq_U2OS_ZNF768_rep2.trim.log
trimmomatic.ChIP_seq_U2OS_ZNF768_rep2.eca5d0f1e9d90118c2dbb58a91b94414.mugqic.done
chmod 755 $COMMAND
trimmomatic_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 20G -c 8 -N 1    | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trimmomatic_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.40b412ebdd74e4613f9afc0388ebbb20.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'merge_trimmomatic_stats.40b412ebdd74e4613f9afc0388ebbb20.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 && \
mkdir -p metrics && \
touch metrics && \

echo -e "Sample\tReadset\tMark Name\tRaw Paired Reads #\tSurviving Paired Reads #\tSurviving Paired Reads %" > metrics/trimReadsetTable.tsv && \
grep ^Input trim/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/RAJI_ZNF768_rep1\tChIP_seq_RAJI_ZNF768_rep1\tZNF768\t\1\t\2/' | \
awk '{OFS="\t"; print $0, $5 / $4 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/RAJI_ZNF768_rep2\tChIP_seq_RAJI_ZNF768_rep2\tZNF768\t\1\t\2/' | \
awk '{OFS="\t"; print $0, $5 / $4 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/U2OS_ZNF768_rep1\tChIP_seq_U2OS_ZNF768_rep1\tZNF768\t\1\t\2/' | \
awk '{OFS="\t"; print $0, $5 / $4 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/U2OS_ZNF768_rep2\tChIP_seq_U2OS_ZNF768_rep2\tZNF768\t\1\t\2/' | \
awk '{OFS="\t"; print $0, $5 / $4 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"\t" '{OFS="\t"; if (NR==1) {if ($3=="Raw Paired Reads #") {paired=1};print "Sample", "Mark Name", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$3=$3*2; $4=$4*2}; sample[$1$2]=$1; markname[$1$2]=$2; raw[$1$2]+=$3; surviving[$1$2]+=$4}}END{for (samplemark in raw){print sample[samplemark], markname[samplemark], raw[samplemark], surviving[samplemark], surviving[samplemark] / raw[samplemark] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/
merge_trimmomatic_stats.40b412ebdd74e4613f9afc0388ebbb20.mugqic.done
chmod 755 $COMMAND
merge_trimmomatic_stats_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=00:20:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: mapping_bwa_mem_sambamba
#-------------------------------------------------------------------------------
STEP=mapping_bwa_mem_sambamba
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: mapping_bwa_mem_sambamba_1_JOB_ID: mapping_bwa_mem_sambamba.ChIP_seq_RAJI_ZNF768_rep1
#-------------------------------------------------------------------------------
JOB_NAME=mapping_bwa_mem_sambamba.ChIP_seq_RAJI_ZNF768_rep1
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/mapping_bwa_mem_sambamba/mapping_bwa_mem_sambamba.ChIP_seq_RAJI_ZNF768_rep1.415989ac5fa061a3646f4388bd3a8c07.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'mapping_bwa_mem_sambamba.ChIP_seq_RAJI_ZNF768_rep1.415989ac5fa061a3646f4388bd3a8c07.mugqic.done' > $COMMAND
module purge && \
module load mugqic/bwa-mem2/2.2.1 mugqic/sambamba/0.8.1 && \
mkdir -p alignment/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1 && \
touch alignment/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1 && \
bwa-mem2 mem -K 100000000 -v 3 -t 16 -Y \
   \
  -R '@RG\tID:ChIP_seq_RAJI_ZNF768_rep1\tSM:RAJI_ZNF768_rep1\tLB:RAJI_ZNF768_rep1\tCN:McGill_Genome_Centre\tPL:Illumina' \
  /cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/bwa-mem2_index/Homo_sapiens.GRCh38.fa \
  trim/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1.trim.pair1.fastq.gz \
  trim/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1.trim.pair2.fastq.gz | \
sambamba view -S -f bam \
  /dev/stdin \
    | \
sambamba sort  \
  /dev/stdin \
  --tmpdir ${SLURM_TMPDIR} \
  --out alignment/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1/ChIP_seq_RAJI_ZNF768_rep1.sorted.bam && \
sambamba index  \
  alignment/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1/ChIP_seq_RAJI_ZNF768_rep1.sorted.bam \
  alignment/RAJI_ZNF768_rep1/ZNF768/ChIP_seq_RAJI_ZNF768_rep1/ChIP_seq_RAJI_ZNF768_rep1.sorted.bam.bai
mapping_bwa_mem_sambamba.ChIP_seq_RAJI_ZNF768_rep1.415989ac5fa061a3646f4388bd3a8c07.mugqic.done
chmod 755 $COMMAND
mapping_bwa_mem_sambamba_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 3900M -c 16 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$mapping_bwa_mem_sambamba_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$mapping_bwa_mem_sambamba_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: mapping_bwa_mem_sambamba_2_JOB_ID: mapping_bwa_mem_sambamba.ChIP_seq_RAJI_ZNF768_rep2
#-------------------------------------------------------------------------------
JOB_NAME=mapping_bwa_mem_sambamba.ChIP_seq_RAJI_ZNF768_rep2
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/mapping_bwa_mem_sambamba/mapping_bwa_mem_sambamba.ChIP_seq_RAJI_ZNF768_rep2.512cda78813640283e12a19f5a67ad45.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'mapping_bwa_mem_sambamba.ChIP_seq_RAJI_ZNF768_rep2.512cda78813640283e12a19f5a67ad45.mugqic.done' > $COMMAND
module purge && \
module load mugqic/bwa-mem2/2.2.1 mugqic/sambamba/0.8.1 && \
mkdir -p alignment/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2 && \
touch alignment/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2 && \
bwa-mem2 mem -K 100000000 -v 3 -t 16 -Y \
   \
  -R '@RG\tID:ChIP_seq_RAJI_ZNF768_rep2\tSM:RAJI_ZNF768_rep2\tLB:RAJI_ZNF768_rep2\tCN:McGill_Genome_Centre\tPL:Illumina' \
  /cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/bwa-mem2_index/Homo_sapiens.GRCh38.fa \
  trim/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2.trim.pair1.fastq.gz \
  trim/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2.trim.pair2.fastq.gz | \
sambamba view -S -f bam \
  /dev/stdin \
    | \
sambamba sort  \
  /dev/stdin \
  --tmpdir ${SLURM_TMPDIR} \
  --out alignment/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2/ChIP_seq_RAJI_ZNF768_rep2.sorted.bam && \
sambamba index  \
  alignment/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2/ChIP_seq_RAJI_ZNF768_rep2.sorted.bam \
  alignment/RAJI_ZNF768_rep2/ZNF768/ChIP_seq_RAJI_ZNF768_rep2/ChIP_seq_RAJI_ZNF768_rep2.sorted.bam.bai
mapping_bwa_mem_sambamba.ChIP_seq_RAJI_ZNF768_rep2.512cda78813640283e12a19f5a67ad45.mugqic.done
chmod 755 $COMMAND
mapping_bwa_mem_sambamba_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 3900M -c 16 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$mapping_bwa_mem_sambamba_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$mapping_bwa_mem_sambamba_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: mapping_bwa_mem_sambamba_3_JOB_ID: mapping_bwa_mem_sambamba.ChIP_seq_U2OS_ZNF768_rep1
#-------------------------------------------------------------------------------
JOB_NAME=mapping_bwa_mem_sambamba.ChIP_seq_U2OS_ZNF768_rep1
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/mapping_bwa_mem_sambamba/mapping_bwa_mem_sambamba.ChIP_seq_U2OS_ZNF768_rep1.4777ffcbcf125f143183bccb25a2dcad.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'mapping_bwa_mem_sambamba.ChIP_seq_U2OS_ZNF768_rep1.4777ffcbcf125f143183bccb25a2dcad.mugqic.done' > $COMMAND
module purge && \
module load mugqic/bwa-mem2/2.2.1 mugqic/sambamba/0.8.1 && \
mkdir -p alignment/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1 && \
touch alignment/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1 && \
bwa-mem2 mem -K 100000000 -v 3 -t 16 -Y \
   \
  -R '@RG\tID:ChIP_seq_U2OS_ZNF768_rep1\tSM:U2OS_ZNF768_rep1\tLB:U2OS_ZNF768_rep1\tCN:McGill_Genome_Centre\tPL:Illumina' \
  /cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/bwa-mem2_index/Homo_sapiens.GRCh38.fa \
  trim/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1.trim.pair1.fastq.gz \
  trim/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1.trim.pair2.fastq.gz | \
sambamba view -S -f bam \
  /dev/stdin \
    | \
sambamba sort  \
  /dev/stdin \
  --tmpdir ${SLURM_TMPDIR} \
  --out alignment/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1/ChIP_seq_U2OS_ZNF768_rep1.sorted.bam && \
sambamba index  \
  alignment/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1/ChIP_seq_U2OS_ZNF768_rep1.sorted.bam \
  alignment/U2OS_ZNF768_rep1/ZNF768/ChIP_seq_U2OS_ZNF768_rep1/ChIP_seq_U2OS_ZNF768_rep1.sorted.bam.bai
mapping_bwa_mem_sambamba.ChIP_seq_U2OS_ZNF768_rep1.4777ffcbcf125f143183bccb25a2dcad.mugqic.done
chmod 755 $COMMAND
mapping_bwa_mem_sambamba_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 3900M -c 16 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$mapping_bwa_mem_sambamba_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$mapping_bwa_mem_sambamba_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: mapping_bwa_mem_sambamba_4_JOB_ID: mapping_bwa_mem_sambamba.ChIP_seq_U2OS_ZNF768_rep2
#-------------------------------------------------------------------------------
JOB_NAME=mapping_bwa_mem_sambamba.ChIP_seq_U2OS_ZNF768_rep2
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/mapping_bwa_mem_sambamba/mapping_bwa_mem_sambamba.ChIP_seq_U2OS_ZNF768_rep2.5bd38e4bf2ccacc2fdd64b179ef4d944.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'mapping_bwa_mem_sambamba.ChIP_seq_U2OS_ZNF768_rep2.5bd38e4bf2ccacc2fdd64b179ef4d944.mugqic.done' > $COMMAND
module purge && \
module load mugqic/bwa-mem2/2.2.1 mugqic/sambamba/0.8.1 && \
mkdir -p alignment/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2 && \
touch alignment/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2 && \
bwa-mem2 mem -K 100000000 -v 3 -t 16 -Y \
   \
  -R '@RG\tID:ChIP_seq_U2OS_ZNF768_rep2\tSM:U2OS_ZNF768_rep2\tLB:U2OS_ZNF768_rep2\tCN:McGill_Genome_Centre\tPL:Illumina' \
  /cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/bwa-mem2_index/Homo_sapiens.GRCh38.fa \
  trim/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2.trim.pair1.fastq.gz \
  trim/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2.trim.pair2.fastq.gz | \
sambamba view -S -f bam \
  /dev/stdin \
    | \
sambamba sort  \
  /dev/stdin \
  --tmpdir ${SLURM_TMPDIR} \
  --out alignment/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2/ChIP_seq_U2OS_ZNF768_rep2.sorted.bam && \
sambamba index  \
  alignment/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2/ChIP_seq_U2OS_ZNF768_rep2.sorted.bam \
  alignment/U2OS_ZNF768_rep2/ZNF768/ChIP_seq_U2OS_ZNF768_rep2/ChIP_seq_U2OS_ZNF768_rep2.sorted.bam.bai
mapping_bwa_mem_sambamba.ChIP_seq_U2OS_ZNF768_rep2.5bd38e4bf2ccacc2fdd64b179ef4d944.mugqic.done
chmod 755 $COMMAND
mapping_bwa_mem_sambamba_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 3900M -c 16 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$mapping_bwa_mem_sambamba_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$mapping_bwa_mem_sambamba_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: sambamba_merge_bam_files
#-------------------------------------------------------------------------------
STEP=sambamba_merge_bam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: sambamba_merge_bam_files_1_JOB_ID: symlink_readset_sample_bam.RAJI_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.RAJI_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$mapping_bwa_mem_sambamba_1_JOB_ID
JOB_DONE=job_output/sambamba_merge_bam_files/symlink_readset_sample_bam.RAJI_ZNF768_rep1.ZNF768.605b6caf53a6ceeb97db857cd5e5254b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'symlink_readset_sample_bam.RAJI_ZNF768_rep1.ZNF768.605b6caf53a6ceeb97db857cd5e5254b.mugqic.done' > $COMMAND
mkdir -p alignment/RAJI_ZNF768_rep1/ZNF768 && \
touch alignment/RAJI_ZNF768_rep1/ZNF768 && \
ln -s -f \
  ChIP_seq_RAJI_ZNF768_rep1/ChIP_seq_RAJI_ZNF768_rep1.sorted.bam \
  alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.bam && \
ln -s -f \
  ChIP_seq_RAJI_ZNF768_rep1/ChIP_seq_RAJI_ZNF768_rep1.sorted.bam.bai \
  alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.bam.bai
symlink_readset_sample_bam.RAJI_ZNF768_rep1.ZNF768.605b6caf53a6ceeb97db857cd5e5254b.mugqic.done
chmod 755 $COMMAND
sambamba_merge_bam_files_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=00:10:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_merge_bam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_merge_bam_files_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_merge_bam_files_2_JOB_ID: symlink_readset_sample_bam.RAJI_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.RAJI_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$mapping_bwa_mem_sambamba_2_JOB_ID
JOB_DONE=job_output/sambamba_merge_bam_files/symlink_readset_sample_bam.RAJI_ZNF768_rep2.ZNF768.45029ebf48c66244b4c9d4026acd7a2f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'symlink_readset_sample_bam.RAJI_ZNF768_rep2.ZNF768.45029ebf48c66244b4c9d4026acd7a2f.mugqic.done' > $COMMAND
mkdir -p alignment/RAJI_ZNF768_rep2/ZNF768 && \
touch alignment/RAJI_ZNF768_rep2/ZNF768 && \
ln -s -f \
  ChIP_seq_RAJI_ZNF768_rep2/ChIP_seq_RAJI_ZNF768_rep2.sorted.bam \
  alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.bam && \
ln -s -f \
  ChIP_seq_RAJI_ZNF768_rep2/ChIP_seq_RAJI_ZNF768_rep2.sorted.bam.bai \
  alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.bam.bai
symlink_readset_sample_bam.RAJI_ZNF768_rep2.ZNF768.45029ebf48c66244b4c9d4026acd7a2f.mugqic.done
chmod 755 $COMMAND
sambamba_merge_bam_files_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=00:10:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_merge_bam_files_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_merge_bam_files_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_merge_bam_files_3_JOB_ID: symlink_readset_sample_bam.U2OS_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.U2OS_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$mapping_bwa_mem_sambamba_3_JOB_ID
JOB_DONE=job_output/sambamba_merge_bam_files/symlink_readset_sample_bam.U2OS_ZNF768_rep1.ZNF768.9333728cbcc7e89af967cc00740e08f9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'symlink_readset_sample_bam.U2OS_ZNF768_rep1.ZNF768.9333728cbcc7e89af967cc00740e08f9.mugqic.done' > $COMMAND
mkdir -p alignment/U2OS_ZNF768_rep1/ZNF768 && \
touch alignment/U2OS_ZNF768_rep1/ZNF768 && \
ln -s -f \
  ChIP_seq_U2OS_ZNF768_rep1/ChIP_seq_U2OS_ZNF768_rep1.sorted.bam \
  alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.bam && \
ln -s -f \
  ChIP_seq_U2OS_ZNF768_rep1/ChIP_seq_U2OS_ZNF768_rep1.sorted.bam.bai \
  alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.bam.bai
symlink_readset_sample_bam.U2OS_ZNF768_rep1.ZNF768.9333728cbcc7e89af967cc00740e08f9.mugqic.done
chmod 755 $COMMAND
sambamba_merge_bam_files_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=00:10:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_merge_bam_files_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_merge_bam_files_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_merge_bam_files_4_JOB_ID: symlink_readset_sample_bam.U2OS_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.U2OS_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$mapping_bwa_mem_sambamba_4_JOB_ID
JOB_DONE=job_output/sambamba_merge_bam_files/symlink_readset_sample_bam.U2OS_ZNF768_rep2.ZNF768.d7e68ba9aae95695368d564f54a84fbd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'symlink_readset_sample_bam.U2OS_ZNF768_rep2.ZNF768.d7e68ba9aae95695368d564f54a84fbd.mugqic.done' > $COMMAND
mkdir -p alignment/U2OS_ZNF768_rep2/ZNF768 && \
touch alignment/U2OS_ZNF768_rep2/ZNF768 && \
ln -s -f \
  ChIP_seq_U2OS_ZNF768_rep2/ChIP_seq_U2OS_ZNF768_rep2.sorted.bam \
  alignment/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.bam && \
ln -s -f \
  ChIP_seq_U2OS_ZNF768_rep2/ChIP_seq_U2OS_ZNF768_rep2.sorted.bam.bai \
  alignment/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.bam.bai
symlink_readset_sample_bam.U2OS_ZNF768_rep2.ZNF768.d7e68ba9aae95695368d564f54a84fbd.mugqic.done
chmod 755 $COMMAND
sambamba_merge_bam_files_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=00:10:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_merge_bam_files_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_merge_bam_files_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: sambamba_mark_duplicates
#-------------------------------------------------------------------------------
STEP=sambamba_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: sambamba_mark_duplicates_1_JOB_ID: sambamba_mark_duplicates.RAJI_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_mark_duplicates.RAJI_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$sambamba_merge_bam_files_1_JOB_ID
JOB_DONE=job_output/sambamba_mark_duplicates/sambamba_mark_duplicates.RAJI_ZNF768_rep1.ZNF768.f126f530235c58187e143b562ffc2676.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_mark_duplicates.RAJI_ZNF768_rep1.ZNF768.f126f530235c58187e143b562ffc2676.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p alignment/RAJI_ZNF768_rep1/ZNF768 && \
touch alignment/RAJI_ZNF768_rep1/ZNF768 && \
sambamba markdup -t 4 --sort-buffer-size=8192 --io-buffer-size=1024 \
  alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.bam \
  --tmpdir ${SLURM_TMPDIR} \
  alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.bam
sambamba_mark_duplicates.RAJI_ZNF768_rep1.ZNF768.f126f530235c58187e143b562ffc2676.mugqic.done
chmod 755 $COMMAND
sambamba_mark_duplicates_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 3900M -c 4 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_mark_duplicates_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_mark_duplicates_2_JOB_ID: sambamba_mark_duplicates.RAJI_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_mark_duplicates.RAJI_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$sambamba_merge_bam_files_2_JOB_ID
JOB_DONE=job_output/sambamba_mark_duplicates/sambamba_mark_duplicates.RAJI_ZNF768_rep2.ZNF768.baffe1126df27fa1ae344e860e2583de.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_mark_duplicates.RAJI_ZNF768_rep2.ZNF768.baffe1126df27fa1ae344e860e2583de.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p alignment/RAJI_ZNF768_rep2/ZNF768 && \
touch alignment/RAJI_ZNF768_rep2/ZNF768 && \
sambamba markdup -t 4 --sort-buffer-size=8192 --io-buffer-size=1024 \
  alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.bam \
  --tmpdir ${SLURM_TMPDIR} \
  alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.bam
sambamba_mark_duplicates.RAJI_ZNF768_rep2.ZNF768.baffe1126df27fa1ae344e860e2583de.mugqic.done
chmod 755 $COMMAND
sambamba_mark_duplicates_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 3900M -c 4 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_mark_duplicates_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_mark_duplicates_3_JOB_ID: sambamba_mark_duplicates.U2OS_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_mark_duplicates.U2OS_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$sambamba_merge_bam_files_3_JOB_ID
JOB_DONE=job_output/sambamba_mark_duplicates/sambamba_mark_duplicates.U2OS_ZNF768_rep1.ZNF768.e5b8a6836c5594700496c33036420706.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_mark_duplicates.U2OS_ZNF768_rep1.ZNF768.e5b8a6836c5594700496c33036420706.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p alignment/U2OS_ZNF768_rep1/ZNF768 && \
touch alignment/U2OS_ZNF768_rep1/ZNF768 && \
sambamba markdup -t 4 --sort-buffer-size=8192 --io-buffer-size=1024 \
  alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.bam \
  --tmpdir ${SLURM_TMPDIR} \
  alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.bam
sambamba_mark_duplicates.U2OS_ZNF768_rep1.ZNF768.e5b8a6836c5594700496c33036420706.mugqic.done
chmod 755 $COMMAND
sambamba_mark_duplicates_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 3900M -c 4 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_mark_duplicates_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_mark_duplicates_4_JOB_ID: sambamba_mark_duplicates.U2OS_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_mark_duplicates.U2OS_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$sambamba_merge_bam_files_4_JOB_ID
JOB_DONE=job_output/sambamba_mark_duplicates/sambamba_mark_duplicates.U2OS_ZNF768_rep2.ZNF768.d7a34f7f117955551110b670ca3692bd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_mark_duplicates.U2OS_ZNF768_rep2.ZNF768.d7a34f7f117955551110b670ca3692bd.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p alignment/U2OS_ZNF768_rep2/ZNF768 && \
touch alignment/U2OS_ZNF768_rep2/ZNF768 && \
sambamba markdup -t 4 --sort-buffer-size=8192 --io-buffer-size=1024 \
  alignment/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.bam \
  --tmpdir ${SLURM_TMPDIR} \
  alignment/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.bam
sambamba_mark_duplicates.U2OS_ZNF768_rep2.ZNF768.d7a34f7f117955551110b670ca3692bd.mugqic.done
chmod 755 $COMMAND
sambamba_mark_duplicates_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 3900M -c 4 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_mark_duplicates_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: sambamba_view_filter
#-------------------------------------------------------------------------------
STEP=sambamba_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: sambamba_view_filter_1_JOB_ID: sambamba_view_filter.RAJI_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_view_filter.RAJI_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$sambamba_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/sambamba_view_filter/sambamba_view_filter.RAJI_ZNF768_rep1.ZNF768.901859c000f86fc77db911e6acb95074.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_view_filter.RAJI_ZNF768_rep1.ZNF768.901859c000f86fc77db911e6acb95074.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p alignment/RAJI_ZNF768_rep1/ZNF768 && \
touch alignment/RAJI_ZNF768_rep1/ZNF768 && \
sambamba view -t 4 -f bam -F "not unmapped and not failed_quality_control and mapping_quality >= 20" \
  alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.bam \
  -o alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam  && \
sambamba index  \
  alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam \
  alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam.bai
sambamba_view_filter.RAJI_ZNF768_rep1.ZNF768.901859c000f86fc77db911e6acb95074.mugqic.done
chmod 755 $COMMAND
sambamba_view_filter_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 3900M -c 4 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_view_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_view_filter_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_view_filter_2_JOB_ID: sambamba_view_filter.RAJI_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_view_filter.RAJI_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$sambamba_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/sambamba_view_filter/sambamba_view_filter.RAJI_ZNF768_rep2.ZNF768.ddeb069c992754e20a531c8c0acc651b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_view_filter.RAJI_ZNF768_rep2.ZNF768.ddeb069c992754e20a531c8c0acc651b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p alignment/RAJI_ZNF768_rep2/ZNF768 && \
touch alignment/RAJI_ZNF768_rep2/ZNF768 && \
sambamba view -t 4 -f bam -F "not unmapped and not failed_quality_control and mapping_quality >= 20" \
  alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.bam \
  -o alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam  && \
sambamba index  \
  alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam \
  alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam.bai
sambamba_view_filter.RAJI_ZNF768_rep2.ZNF768.ddeb069c992754e20a531c8c0acc651b.mugqic.done
chmod 755 $COMMAND
sambamba_view_filter_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 3900M -c 4 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_view_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_view_filter_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_view_filter_3_JOB_ID: sambamba_view_filter.U2OS_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_view_filter.U2OS_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$sambamba_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/sambamba_view_filter/sambamba_view_filter.U2OS_ZNF768_rep1.ZNF768.eba845fa8a23f70ca8baa88778fec9f3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_view_filter.U2OS_ZNF768_rep1.ZNF768.eba845fa8a23f70ca8baa88778fec9f3.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p alignment/U2OS_ZNF768_rep1/ZNF768 && \
touch alignment/U2OS_ZNF768_rep1/ZNF768 && \
sambamba view -t 4 -f bam -F "not unmapped and not failed_quality_control and mapping_quality >= 20" \
  alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.bam \
  -o alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam  && \
sambamba index  \
  alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam \
  alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam.bai
sambamba_view_filter.U2OS_ZNF768_rep1.ZNF768.eba845fa8a23f70ca8baa88778fec9f3.mugqic.done
chmod 755 $COMMAND
sambamba_view_filter_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 3900M -c 4 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_view_filter_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_view_filter_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_view_filter_4_JOB_ID: sambamba_view_filter.U2OS_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_view_filter.U2OS_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$sambamba_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/sambamba_view_filter/sambamba_view_filter.U2OS_ZNF768_rep2.ZNF768.a3e5bf671097913002ebb64916c2e9a2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_view_filter.U2OS_ZNF768_rep2.ZNF768.a3e5bf671097913002ebb64916c2e9a2.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p alignment/U2OS_ZNF768_rep2/ZNF768 && \
touch alignment/U2OS_ZNF768_rep2/ZNF768 && \
sambamba view -t 4 -f bam -F "not unmapped and not failed_quality_control and mapping_quality >= 20" \
  alignment/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.bam \
  -o alignment/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam  && \
sambamba index  \
  alignment/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam \
  alignment/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam.bai
sambamba_view_filter.U2OS_ZNF768_rep2.ZNF768.a3e5bf671097913002ebb64916c2e9a2.mugqic.done
chmod 755 $COMMAND
sambamba_view_filter_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 3900M -c 4 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_view_filter_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_view_filter_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: metrics
#-------------------------------------------------------------------------------
STEP=metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: metrics_1_JOB_ID: picard_collect_multiple_metrics.RAJI_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=picard_collect_multiple_metrics.RAJI_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$sambamba_view_filter_1_JOB_ID
JOB_DONE=job_output/metrics/picard_collect_multiple_metrics.RAJI_ZNF768_rep1.ZNF768.3da1b2d6c38a61e566e1a5fea8885729.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_collect_multiple_metrics.RAJI_ZNF768_rep1.ZNF768.3da1b2d6c38a61e566e1a5fea8885729.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.26.6 mugqic/R_Bioconductor/4.0.3_3.12 && \
mkdir -p metrics/RAJI_ZNF768_rep1/ZNF768 && \
touch metrics/RAJI_ZNF768_rep1/ZNF768 && \
mkdir -p metrics/multiqc_inputs && \
touch metrics/multiqc_inputs && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx8G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam \
 OUTPUT=metrics/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics \
 MAX_RECORDS_IN_RAM=1000000 && \
ln -s -f \
  ../RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle_metrics \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle_metrics && \
ln -s -f \
  ../RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle.pdf \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle.pdf && \
ln -s -f \
  ../RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.alignment_summary_metrics \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.alignment_summary_metrics && \
ln -s -f \
  ../RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.insert_size_histogram.pdf \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.insert_size_histogram.pdf && \
ln -s -f \
  ../RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.insert_size_metrics \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.insert_size_metrics && \
ln -s -f \
  ../RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle_metrics \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle_metrics && \
ln -s -f \
  ../RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle.pdf \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle.pdf && \
ln -s -f \
  ../RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution_metrics \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution_metrics && \
ln -s -f \
  ../RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution.pdf \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution.pdf
picard_collect_multiple_metrics.RAJI_ZNF768_rep1.ZNF768.3da1b2d6c38a61e566e1a5fea8885729.mugqic.done
chmod 755 $COMMAND
metrics_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=06:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: metrics_2_JOB_ID: metrics_flagstat.RAJI_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=metrics_flagstat.RAJI_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$sambamba_mark_duplicates_1_JOB_ID:$sambamba_view_filter_1_JOB_ID
JOB_DONE=job_output/metrics/metrics_flagstat.RAJI_ZNF768_rep1.ZNF768.6de13d0b0d176bac9e1a0bb2159250ff.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'metrics_flagstat.RAJI_ZNF768_rep1.ZNF768.6de13d0b0d176bac9e1a0bb2159250ff.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p metrics/RAJI_ZNF768_rep1/ZNF768 && \
touch metrics/RAJI_ZNF768_rep1/ZNF768 && \
sambamba flagstat  \
  alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.bam \
  > metrics/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.flagstat && \
sambamba flagstat  \
  alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam \
  > metrics/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.flagstat && \
ln -s -f \
  ../RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.flagstat \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768.sorted.dup.flagstat
metrics_flagstat.RAJI_ZNF768_rep1.ZNF768.6de13d0b0d176bac9e1a0bb2159250ff.mugqic.done
chmod 755 $COMMAND
metrics_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=06:00:00 --mem 8G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: metrics_3_JOB_ID: picard_collect_multiple_metrics.RAJI_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=picard_collect_multiple_metrics.RAJI_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$sambamba_view_filter_2_JOB_ID
JOB_DONE=job_output/metrics/picard_collect_multiple_metrics.RAJI_ZNF768_rep2.ZNF768.f31946ec463a485ef55bac84198b0994.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_collect_multiple_metrics.RAJI_ZNF768_rep2.ZNF768.f31946ec463a485ef55bac84198b0994.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.26.6 mugqic/R_Bioconductor/4.0.3_3.12 && \
mkdir -p metrics/RAJI_ZNF768_rep2/ZNF768 && \
touch metrics/RAJI_ZNF768_rep2/ZNF768 && \
mkdir -p metrics/multiqc_inputs && \
touch metrics/multiqc_inputs && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx8G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam \
 OUTPUT=metrics/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics \
 MAX_RECORDS_IN_RAM=1000000 && \
ln -s -f \
  ../RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle_metrics \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle_metrics && \
ln -s -f \
  ../RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle.pdf \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle.pdf && \
ln -s -f \
  ../RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.alignment_summary_metrics \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.alignment_summary_metrics && \
ln -s -f \
  ../RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.insert_size_histogram.pdf \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.insert_size_histogram.pdf && \
ln -s -f \
  ../RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.insert_size_metrics \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.insert_size_metrics && \
ln -s -f \
  ../RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle_metrics \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle_metrics && \
ln -s -f \
  ../RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle.pdf \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle.pdf && \
ln -s -f \
  ../RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution_metrics \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution_metrics && \
ln -s -f \
  ../RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution.pdf \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution.pdf
picard_collect_multiple_metrics.RAJI_ZNF768_rep2.ZNF768.f31946ec463a485ef55bac84198b0994.mugqic.done
chmod 755 $COMMAND
metrics_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=06:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: metrics_4_JOB_ID: metrics_flagstat.RAJI_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=metrics_flagstat.RAJI_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$sambamba_mark_duplicates_2_JOB_ID:$sambamba_view_filter_2_JOB_ID
JOB_DONE=job_output/metrics/metrics_flagstat.RAJI_ZNF768_rep2.ZNF768.b04b56218d329d13e8a9860247989359.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'metrics_flagstat.RAJI_ZNF768_rep2.ZNF768.b04b56218d329d13e8a9860247989359.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p metrics/RAJI_ZNF768_rep2/ZNF768 && \
touch metrics/RAJI_ZNF768_rep2/ZNF768 && \
sambamba flagstat  \
  alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.bam \
  > metrics/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.flagstat && \
sambamba flagstat  \
  alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam \
  > metrics/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.flagstat && \
ln -s -f \
  ../RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.flagstat \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768.sorted.dup.flagstat
metrics_flagstat.RAJI_ZNF768_rep2.ZNF768.b04b56218d329d13e8a9860247989359.mugqic.done
chmod 755 $COMMAND
metrics_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=06:00:00 --mem 8G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: metrics_5_JOB_ID: picard_collect_multiple_metrics.U2OS_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=picard_collect_multiple_metrics.U2OS_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$sambamba_view_filter_3_JOB_ID
JOB_DONE=job_output/metrics/picard_collect_multiple_metrics.U2OS_ZNF768_rep1.ZNF768.d2e9d3a39fa45eaf1c9fab2ca1ad363b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_collect_multiple_metrics.U2OS_ZNF768_rep1.ZNF768.d2e9d3a39fa45eaf1c9fab2ca1ad363b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.26.6 mugqic/R_Bioconductor/4.0.3_3.12 && \
mkdir -p metrics/U2OS_ZNF768_rep1/ZNF768 && \
touch metrics/U2OS_ZNF768_rep1/ZNF768 && \
mkdir -p metrics/multiqc_inputs && \
touch metrics/multiqc_inputs && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx8G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam \
 OUTPUT=metrics/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics \
 MAX_RECORDS_IN_RAM=1000000 && \
ln -s -f \
  ../U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle_metrics \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle_metrics && \
ln -s -f \
  ../U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle.pdf \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle.pdf && \
ln -s -f \
  ../U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.alignment_summary_metrics \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.alignment_summary_metrics && \
ln -s -f \
  ../U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.insert_size_histogram.pdf \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.insert_size_histogram.pdf && \
ln -s -f \
  ../U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.insert_size_metrics \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.insert_size_metrics && \
ln -s -f \
  ../U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle_metrics \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle_metrics && \
ln -s -f \
  ../U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle.pdf \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle.pdf && \
ln -s -f \
  ../U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution_metrics \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution_metrics && \
ln -s -f \
  ../U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution.pdf \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution.pdf
picard_collect_multiple_metrics.U2OS_ZNF768_rep1.ZNF768.d2e9d3a39fa45eaf1c9fab2ca1ad363b.mugqic.done
chmod 755 $COMMAND
metrics_5_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=06:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_5_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: metrics_6_JOB_ID: metrics_flagstat.U2OS_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=metrics_flagstat.U2OS_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$sambamba_mark_duplicates_3_JOB_ID:$sambamba_view_filter_3_JOB_ID
JOB_DONE=job_output/metrics/metrics_flagstat.U2OS_ZNF768_rep1.ZNF768.b53276df25d000ef4a71411f9647d644.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'metrics_flagstat.U2OS_ZNF768_rep1.ZNF768.b53276df25d000ef4a71411f9647d644.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p metrics/U2OS_ZNF768_rep1/ZNF768 && \
touch metrics/U2OS_ZNF768_rep1/ZNF768 && \
sambamba flagstat  \
  alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.bam \
  > metrics/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.flagstat && \
sambamba flagstat  \
  alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam \
  > metrics/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.flagstat && \
ln -s -f \
  ../U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.flagstat \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768.sorted.dup.flagstat
metrics_flagstat.U2OS_ZNF768_rep1.ZNF768.b53276df25d000ef4a71411f9647d644.mugqic.done
chmod 755 $COMMAND
metrics_6_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=06:00:00 --mem 8G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_6_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: metrics_7_JOB_ID: picard_collect_multiple_metrics.U2OS_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=picard_collect_multiple_metrics.U2OS_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$sambamba_view_filter_4_JOB_ID
JOB_DONE=job_output/metrics/picard_collect_multiple_metrics.U2OS_ZNF768_rep2.ZNF768.f669d58bb981e399a71727ed9742e0fb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_collect_multiple_metrics.U2OS_ZNF768_rep2.ZNF768.f669d58bb981e399a71727ed9742e0fb.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.26.6 mugqic/R_Bioconductor/4.0.3_3.12 && \
mkdir -p metrics/U2OS_ZNF768_rep2/ZNF768 && \
touch metrics/U2OS_ZNF768_rep2/ZNF768 && \
mkdir -p metrics/multiqc_inputs && \
touch metrics/multiqc_inputs && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx8G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam \
 OUTPUT=metrics/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics \
 MAX_RECORDS_IN_RAM=1000000 && \
ln -s -f \
  ../U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle_metrics \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle_metrics && \
ln -s -f \
  ../U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle.pdf \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.base_distribution_by_cycle.pdf && \
ln -s -f \
  ../U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.alignment_summary_metrics \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.alignment_summary_metrics && \
ln -s -f \
  ../U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.insert_size_histogram.pdf \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.insert_size_histogram.pdf && \
ln -s -f \
  ../U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.insert_size_metrics \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.insert_size_metrics && \
ln -s -f \
  ../U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle_metrics \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle_metrics && \
ln -s -f \
  ../U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle.pdf \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_by_cycle.pdf && \
ln -s -f \
  ../U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution_metrics \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution_metrics && \
ln -s -f \
  ../U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution.pdf \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.all.metrics.quality_distribution.pdf
picard_collect_multiple_metrics.U2OS_ZNF768_rep2.ZNF768.f669d58bb981e399a71727ed9742e0fb.mugqic.done
chmod 755 $COMMAND
metrics_7_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=06:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_7_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: metrics_8_JOB_ID: metrics_flagstat.U2OS_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=metrics_flagstat.U2OS_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$sambamba_mark_duplicates_4_JOB_ID:$sambamba_view_filter_4_JOB_ID
JOB_DONE=job_output/metrics/metrics_flagstat.U2OS_ZNF768_rep2.ZNF768.f2da477079a418b38f2c4ac3cfa8ba69.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'metrics_flagstat.U2OS_ZNF768_rep2.ZNF768.f2da477079a418b38f2c4ac3cfa8ba69.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p metrics/U2OS_ZNF768_rep2/ZNF768 && \
touch metrics/U2OS_ZNF768_rep2/ZNF768 && \
sambamba flagstat  \
  alignment/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.bam \
  > metrics/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.flagstat && \
sambamba flagstat  \
  alignment/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam \
  > metrics/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.flagstat && \
ln -s -f \
  ../U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.flagstat \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768.sorted.dup.flagstat
metrics_flagstat.U2OS_ZNF768_rep2.ZNF768.f2da477079a418b38f2c4ac3cfa8ba69.mugqic.done
chmod 755 $COMMAND
metrics_8_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=06:00:00 --mem 8G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_8_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: metrics_9_JOB_ID: metrics_report
#-------------------------------------------------------------------------------
JOB_NAME=metrics_report
JOB_DEPENDENCIES=$sambamba_view_filter_1_JOB_ID:$sambamba_view_filter_2_JOB_ID:$sambamba_view_filter_3_JOB_ID:$sambamba_view_filter_4_JOB_ID:$metrics_2_JOB_ID:$metrics_4_JOB_ID:$metrics_6_JOB_ID:$metrics_8_JOB_ID
JOB_DONE=job_output/metrics/metrics_report.7fe26ddbc9028bad3d0d908918bafdea.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'metrics_report.7fe26ddbc9028bad3d0d908918bafdea.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 mugqic/samtools/1.14 && \
mkdir -p metrics
cp /dev/null metrics/SampleMetrics.tsv && \
declare -A samples_associative_array=(["RAJI_ZNF768_rep1"]="ZNF768" ["RAJI_ZNF768_rep2"]="ZNF768" ["U2OS_ZNF768_rep1"]="ZNF768" ["U2OS_ZNF768_rep2"]="ZNF768") && \
for sample in ${!samples_associative_array[@]}
do
  for mark_name in ${samples_associative_array[$sample]}
  do
    raw_flagstat_file=metrics/$sample/$mark_name/$sample.$mark_name.sorted.dup.flagstat
    filtered_flagstat_file=metrics/$sample/$mark_name/$sample.$mark_name.sorted.dup.filtered.flagstat
    bam_file=alignment/$sample/$mark_name/$sample.$mark_name.sorted.dup.filtered.bam
    raw_supplementarysecondary_reads=`bc <<< $(grep "secondary" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* secondary.*//')+$(grep "supplementary" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* supplementary.*//')`
    mapped_reads=`bc <<< $(grep "mapped (" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* mapped (.*)//')-$raw_supplementarysecondary_reads`
    filtered_supplementarysecondary_reads=`bc <<< $(grep "secondary" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* secondary.*//')+$(grep "supplementary" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* supplementary.*//')`
    filtered_reads=`bc <<< $(grep "in total" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$filtered_supplementarysecondary_reads`
    filtered_mapped_reads=`bc <<< $(grep "mapped (" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* mapped (.*)//')-$filtered_supplementarysecondary_reads`
    filtered_mapped_rate=`echo "scale=4; 100*$filtered_mapped_reads/$filtered_reads" | bc -l`
    filtered_dup_reads=`grep "duplicates" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* duplicates$//'`
    filtered_dup_rate=`echo "scale=4; 100*$filtered_dup_reads/$filtered_mapped_reads" | bc -l`
    filtered_dedup_reads=`echo "$filtered_mapped_reads-$filtered_dup_reads" | bc -l`
    if [[ -s metrics/trimSampleTable.tsv ]]
      then
        raw_reads=$(grep -P "${sample}\t${mark_name}" metrics/trimSampleTable.tsv | cut -f 3)
        raw_trimmed_reads=`bc <<< $(grep "in total" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$raw_supplementarysecondary_reads`
        mapped_reads_rate=`echo "scale=4; 100*$mapped_reads/$raw_trimmed_reads" | bc -l`
        raw_trimmed_rate=`echo "scale=4; 100*$raw_trimmed_reads/$raw_reads" | bc -l`
        filtered_rate=`echo "scale=4; 100*$filtered_reads/$raw_trimmed_reads" | bc -l`
      else
        raw_reads=`bc <<< $(grep "in total" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$raw_supplementarysecondary_reads`
        raw_trimmed_reads="NULL"
        mapped_reads_rate=`echo "scale=4; 100*$mapped_reads/$raw_reads" | bc -l`
        raw_trimmed_rate="NULL"
        filtered_rate=`echo "scale=4; 100*$filtered_reads/$raw_reads" | bc -l`
    fi
    filtered_mito_reads=$(sambamba view -F "not duplicate" -c $bam_file chrM)
    filtered_mito_rate=$(echo "scale=4; 100*$filtered_mito_reads/$filtered_mapped_reads" | bc -l)
    if [[ $mark_name != "Input" ]]
        then
          chip_bed_file=peak_call/$sample/$mark_name/$sample.${mark_name}_peaks.*Peak.bed
          nmb_peaks=$(wc -l $chip_bed_file | cut -f 1 -d " ")
          reads_under_peaks=$(samtools view -c -L $chip_bed_file $bam_file)
          frip=$(echo "scale=4; $reads_under_peaks/$filtered_mapped_reads" | bc -l)
        else
          nmb_peaks="NA"
          reads_under_peaks="NA"
          frip="NA"
    fi
    echo -e "$sample\t$mark_name\t$raw_reads\t$raw_trimmed_reads\t$raw_trimmed_rate\t$mapped_reads\t$mapped_reads_rate\t$filtered_reads\t$filtered_rate\t$filtered_mapped_reads\t$filtered_mapped_rate\t$filtered_dup_reads\t$filtered_dup_rate\t$filtered_dedup_reads\t$filtered_mito_reads\t$filtered_mito_rate\t$nmb_peaks\t$reads_under_peaks\t$frip" >> metrics/SampleMetrics.tsv
  done
done && \
sed -i -e "1 i\Sample\tMark Name\tRaw Reads #\tRemaining Reads after Trimming #\tRemaining Reads after Trimming %\tAligned Trimmed Reads #\tAligned Trimmed Reads %\tRemaining Reads after Filtering #\tRemaining Reads after Filtering %\tAligned Filtered Reads #\tAligned Filtered Reads %\tDuplicate Reads #\tDuplicate Reads %\tFinal Aligned Reads # without Duplicates\tMitochondrial Reads #\tMitochondrial Reads %\tNumber of Peaks\tReads under Peaks\tFRIP" metrics/SampleMetrics.tsv && \
mkdir -p report && \
cp metrics/SampleMetrics.tsv report/SampleMetrics.tsv
metrics_report.7fe26ddbc9028bad3d0d908918bafdea.mugqic.done
chmod 755 $COMMAND
metrics_9_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=03:00:00 --mem 16G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_9_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: homer_make_tag_directory
#-------------------------------------------------------------------------------
STEP=homer_make_tag_directory
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.RAJI_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.RAJI_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$sambamba_view_filter_1_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.RAJI_ZNF768_rep1.ZNF768.bcbf17c0ba95658da91da3771f0f8b26.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_tag_directory.RAJI_ZNF768_rep1.ZNF768.bcbf17c0ba95658da91da3771f0f8b26.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.11 mugqic/samtools/1.14 && \
mkdir -p metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768 && \
touch metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768 && \
makeTagDirectory \
  tags/RAJI_ZNF768_rep1/RAJI_ZNF768_rep1.ZNF768 \
  alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam \
  -genome hg38 \
  -checkGC && \
ln -s -f \
  ../../../tags/RAJI_ZNF768_rep1/RAJI_ZNF768_rep1.ZNF768/tagInfo.txt \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768/tagInfo.txt && \
ln -s -f \
  ../../../tags/RAJI_ZNF768_rep1/RAJI_ZNF768_rep1.ZNF768/tagGCcontent.txt \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768/tagGCcontent.txt && \
ln -s -f \
  ../../../tags/RAJI_ZNF768_rep1/RAJI_ZNF768_rep1.ZNF768/genomeGCcontent.txt \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768/genomeGCcontent.txt && \
ln -s -f \
  ../../../tags/RAJI_ZNF768_rep1/RAJI_ZNF768_rep1.ZNF768/tagLengthDistribution.txt \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768/tagLengthDistribution.txt
homer_make_tag_directory.RAJI_ZNF768_rep1.ZNF768.bcbf17c0ba95658da91da3771f0f8b26.mugqic.done
chmod 755 $COMMAND
homer_make_tag_directory_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 8G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.RAJI_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.RAJI_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$sambamba_view_filter_2_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.RAJI_ZNF768_rep2.ZNF768.1c480bfa6ee747ff572055f46b47309f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_tag_directory.RAJI_ZNF768_rep2.ZNF768.1c480bfa6ee747ff572055f46b47309f.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.11 mugqic/samtools/1.14 && \
mkdir -p metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768 && \
touch metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768 && \
makeTagDirectory \
  tags/RAJI_ZNF768_rep2/RAJI_ZNF768_rep2.ZNF768 \
  alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam \
  -genome hg38 \
  -checkGC && \
ln -s -f \
  ../../../tags/RAJI_ZNF768_rep2/RAJI_ZNF768_rep2.ZNF768/tagInfo.txt \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768/tagInfo.txt && \
ln -s -f \
  ../../../tags/RAJI_ZNF768_rep2/RAJI_ZNF768_rep2.ZNF768/tagGCcontent.txt \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768/tagGCcontent.txt && \
ln -s -f \
  ../../../tags/RAJI_ZNF768_rep2/RAJI_ZNF768_rep2.ZNF768/genomeGCcontent.txt \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768/genomeGCcontent.txt && \
ln -s -f \
  ../../../tags/RAJI_ZNF768_rep2/RAJI_ZNF768_rep2.ZNF768/tagLengthDistribution.txt \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768/tagLengthDistribution.txt
homer_make_tag_directory.RAJI_ZNF768_rep2.ZNF768.1c480bfa6ee747ff572055f46b47309f.mugqic.done
chmod 755 $COMMAND
homer_make_tag_directory_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 8G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.U2OS_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.U2OS_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$sambamba_view_filter_3_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.U2OS_ZNF768_rep1.ZNF768.382f29deb90133266af4589ff3af49a6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_tag_directory.U2OS_ZNF768_rep1.ZNF768.382f29deb90133266af4589ff3af49a6.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.11 mugqic/samtools/1.14 && \
mkdir -p metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768 && \
touch metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768 && \
makeTagDirectory \
  tags/U2OS_ZNF768_rep1/U2OS_ZNF768_rep1.ZNF768 \
  alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam \
  -genome hg38 \
  -checkGC && \
ln -s -f \
  ../../../tags/U2OS_ZNF768_rep1/U2OS_ZNF768_rep1.ZNF768/tagInfo.txt \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768/tagInfo.txt && \
ln -s -f \
  ../../../tags/U2OS_ZNF768_rep1/U2OS_ZNF768_rep1.ZNF768/tagGCcontent.txt \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768/tagGCcontent.txt && \
ln -s -f \
  ../../../tags/U2OS_ZNF768_rep1/U2OS_ZNF768_rep1.ZNF768/genomeGCcontent.txt \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768/genomeGCcontent.txt && \
ln -s -f \
  ../../../tags/U2OS_ZNF768_rep1/U2OS_ZNF768_rep1.ZNF768/tagLengthDistribution.txt \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768/tagLengthDistribution.txt
homer_make_tag_directory.U2OS_ZNF768_rep1.ZNF768.382f29deb90133266af4589ff3af49a6.mugqic.done
chmod 755 $COMMAND
homer_make_tag_directory_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 8G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.U2OS_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.U2OS_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$sambamba_view_filter_4_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.U2OS_ZNF768_rep2.ZNF768.089bb6a0b76d1cf7821f99f0d2d5d3ce.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_tag_directory.U2OS_ZNF768_rep2.ZNF768.089bb6a0b76d1cf7821f99f0d2d5d3ce.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.11 mugqic/samtools/1.14 && \
mkdir -p metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768 && \
touch metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768 && \
makeTagDirectory \
  tags/U2OS_ZNF768_rep2/U2OS_ZNF768_rep2.ZNF768 \
  alignment/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam \
  -genome hg38 \
  -checkGC && \
ln -s -f \
  ../../../tags/U2OS_ZNF768_rep2/U2OS_ZNF768_rep2.ZNF768/tagInfo.txt \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768/tagInfo.txt && \
ln -s -f \
  ../../../tags/U2OS_ZNF768_rep2/U2OS_ZNF768_rep2.ZNF768/tagGCcontent.txt \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768/tagGCcontent.txt && \
ln -s -f \
  ../../../tags/U2OS_ZNF768_rep2/U2OS_ZNF768_rep2.ZNF768/genomeGCcontent.txt \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768/genomeGCcontent.txt && \
ln -s -f \
  ../../../tags/U2OS_ZNF768_rep2/U2OS_ZNF768_rep2.ZNF768/tagLengthDistribution.txt \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768/tagLengthDistribution.txt
homer_make_tag_directory.U2OS_ZNF768_rep2.ZNF768.089bb6a0b76d1cf7821f99f0d2d5d3ce.mugqic.done
chmod 755 $COMMAND
homer_make_tag_directory_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 8G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_tag_directory_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_tag_directory_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: qc_metrics
#-------------------------------------------------------------------------------
STEP=qc_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: qc_metrics_1_JOB_ID: qc_plots_R
#-------------------------------------------------------------------------------
JOB_NAME=qc_plots_R
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID:$homer_make_tag_directory_2_JOB_ID:$homer_make_tag_directory_3_JOB_ID:$homer_make_tag_directory_4_JOB_ID
JOB_DONE=job_output/qc_metrics/qc_plots_R.379cd63ce233073def46a8a5758cbee8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'qc_plots_R.379cd63ce233073def46a8a5758cbee8.mugqic.done' > $COMMAND
module purge && \
module load mugqic/mugqic_tools/2.12.4 mugqic/R_Bioconductor/4.0.3_3.12 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/chipseq_ZNF768_GSE111879/readset_chipseq_ZNF768_GSE111879_20250221.txt \
  /lustre06/project/6001942/chris11/20250221_ZNF768_project/output/chip-pipeline_ZNF768_GSE111879-GRCh38_PE && \
declare -A samples_associative_array=(["RAJI_ZNF768_rep1"]="ZNF768" ["RAJI_ZNF768_rep2"]="ZNF768" ["U2OS_ZNF768_rep1"]="ZNF768" ["U2OS_ZNF768_rep2"]="ZNF768") && \
for sample in ${!samples_associative_array[@]}
do
  for mark_name in ${samples_associative_array[$sample]}
  do
    cp --parents graphs/${sample}.${mark_name}_QC_Metrics.ps report/
    convert -rotate 90 graphs/${sample}.${mark_name}_QC_Metrics.ps report/graphs/${sample}.${mark_name}_QC_Metrics.png
  done
done
qc_plots_R.379cd63ce233073def46a8a5758cbee8.mugqic.done
chmod 755 $COMMAND
qc_metrics_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$qc_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$qc_metrics_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: homer_make_ucsc_file
#-------------------------------------------------------------------------------
STEP=homer_make_ucsc_file
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_1_JOB_ID: homer_make_ucsc_file.RAJI_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.RAJI_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.RAJI_ZNF768_rep1.ZNF768.cd27abf12d7b5d975aa96724b302160b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file.RAJI_ZNF768_rep1.ZNF768.cd27abf12d7b5d975aa96724b302160b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.11 && \
mkdir -p tracks/RAJI_ZNF768_rep1/ZNF768 && \
touch tracks/RAJI_ZNF768_rep1/ZNF768 && \
makeUCSCfile \
  tags/RAJI_ZNF768_rep1/RAJI_ZNF768_rep1.ZNF768 \
  > tracks/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.ucsc.bedGraph && \
gzip -c tracks/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.ucsc.bedGraph \
  > tracks/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.ucsc.bedGraph.gz
homer_make_ucsc_file.RAJI_ZNF768_rep1.ZNF768.cd27abf12d7b5d975aa96724b302160b.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 32G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_2_JOB_ID: homer_make_ucsc_file_bigWig.RAJI_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.RAJI_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$homer_make_ucsc_file_1_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.RAJI_ZNF768_rep1.ZNF768.ca47069d97a257f7080e0c242bed0e31.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file_bigWig.RAJI_ZNF768_rep1.ZNF768.ca47069d97a257f7080e0c242bed0e31.mugqic.done' > $COMMAND
module purge && \
module load mugqic/ucsc/v346 && \
mkdir -p tracks/RAJI_ZNF768_rep1/ZNF768/bigWig && \
touch tracks/RAJI_ZNF768_rep1/ZNF768/bigWig && \
(cat tracks/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.ucsc.bedGraph | head -n 1 > tracks/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.ucsc.bedGraph.head.tmp ; ec=$?; if [ "$ec" -eq 141 ]; then exit 0; else exit "$ec"; fi) && \
cat tracks/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -v "GL\|lambda\|pUC19\|KI\|\KN\|random"  | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.ucsc.bedGraph.body.tmp && \
cat tracks/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.ucsc.bedGraph.head.tmp tracks/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.ucsc.bedGraph.body.tmp > tracks/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.ucsc.bedGraph.sorted && \
rm tracks/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.ucsc.bedGraph.head.tmp tracks/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/RAJI_ZNF768_rep1/ZNF768/bigWig/RAJI_ZNF768_rep1.ZNF768.bw
homer_make_ucsc_file_bigWig.RAJI_ZNF768_rep1.ZNF768.ca47069d97a257f7080e0c242bed0e31.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_3_JOB_ID: homer_make_ucsc_file.RAJI_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.RAJI_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$homer_make_tag_directory_2_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.RAJI_ZNF768_rep2.ZNF768.020b5e9cdbe6e094ec2e0ba6c126bfe4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file.RAJI_ZNF768_rep2.ZNF768.020b5e9cdbe6e094ec2e0ba6c126bfe4.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.11 && \
mkdir -p tracks/RAJI_ZNF768_rep2/ZNF768 && \
touch tracks/RAJI_ZNF768_rep2/ZNF768 && \
makeUCSCfile \
  tags/RAJI_ZNF768_rep2/RAJI_ZNF768_rep2.ZNF768 \
  > tracks/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.ucsc.bedGraph && \
gzip -c tracks/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.ucsc.bedGraph \
  > tracks/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.ucsc.bedGraph.gz
homer_make_ucsc_file.RAJI_ZNF768_rep2.ZNF768.020b5e9cdbe6e094ec2e0ba6c126bfe4.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 32G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_4_JOB_ID: homer_make_ucsc_file_bigWig.RAJI_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.RAJI_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$homer_make_ucsc_file_3_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.RAJI_ZNF768_rep2.ZNF768.348e5800542aaa3c7ea051a0f582e0c4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file_bigWig.RAJI_ZNF768_rep2.ZNF768.348e5800542aaa3c7ea051a0f582e0c4.mugqic.done' > $COMMAND
module purge && \
module load mugqic/ucsc/v346 && \
mkdir -p tracks/RAJI_ZNF768_rep2/ZNF768/bigWig && \
touch tracks/RAJI_ZNF768_rep2/ZNF768/bigWig && \
(cat tracks/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.ucsc.bedGraph | head -n 1 > tracks/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.ucsc.bedGraph.head.tmp ; ec=$?; if [ "$ec" -eq 141 ]; then exit 0; else exit "$ec"; fi) && \
cat tracks/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -v "GL\|lambda\|pUC19\|KI\|\KN\|random"  | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.ucsc.bedGraph.body.tmp && \
cat tracks/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.ucsc.bedGraph.head.tmp tracks/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.ucsc.bedGraph.body.tmp > tracks/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.ucsc.bedGraph.sorted && \
rm tracks/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.ucsc.bedGraph.head.tmp tracks/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/RAJI_ZNF768_rep2/ZNF768/bigWig/RAJI_ZNF768_rep2.ZNF768.bw
homer_make_ucsc_file_bigWig.RAJI_ZNF768_rep2.ZNF768.348e5800542aaa3c7ea051a0f582e0c4.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_5_JOB_ID: homer_make_ucsc_file.U2OS_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.U2OS_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$homer_make_tag_directory_3_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.U2OS_ZNF768_rep1.ZNF768.02d68fab2b479cd3a1b9c3fc8ecdc481.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file.U2OS_ZNF768_rep1.ZNF768.02d68fab2b479cd3a1b9c3fc8ecdc481.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.11 && \
mkdir -p tracks/U2OS_ZNF768_rep1/ZNF768 && \
touch tracks/U2OS_ZNF768_rep1/ZNF768 && \
makeUCSCfile \
  tags/U2OS_ZNF768_rep1/U2OS_ZNF768_rep1.ZNF768 \
  > tracks/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.ucsc.bedGraph && \
gzip -c tracks/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.ucsc.bedGraph \
  > tracks/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.ucsc.bedGraph.gz
homer_make_ucsc_file.U2OS_ZNF768_rep1.ZNF768.02d68fab2b479cd3a1b9c3fc8ecdc481.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_5_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 32G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_5_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_6_JOB_ID: homer_make_ucsc_file_bigWig.U2OS_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.U2OS_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$homer_make_ucsc_file_5_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.U2OS_ZNF768_rep1.ZNF768.fb6fd9c14ec0650daad095d05a8fb9cd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file_bigWig.U2OS_ZNF768_rep1.ZNF768.fb6fd9c14ec0650daad095d05a8fb9cd.mugqic.done' > $COMMAND
module purge && \
module load mugqic/ucsc/v346 && \
mkdir -p tracks/U2OS_ZNF768_rep1/ZNF768/bigWig && \
touch tracks/U2OS_ZNF768_rep1/ZNF768/bigWig && \
(cat tracks/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.ucsc.bedGraph | head -n 1 > tracks/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.ucsc.bedGraph.head.tmp ; ec=$?; if [ "$ec" -eq 141 ]; then exit 0; else exit "$ec"; fi) && \
cat tracks/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -v "GL\|lambda\|pUC19\|KI\|\KN\|random"  | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.ucsc.bedGraph.body.tmp && \
cat tracks/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.ucsc.bedGraph.head.tmp tracks/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.ucsc.bedGraph.body.tmp > tracks/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.ucsc.bedGraph.sorted && \
rm tracks/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.ucsc.bedGraph.head.tmp tracks/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/U2OS_ZNF768_rep1/ZNF768/bigWig/U2OS_ZNF768_rep1.ZNF768.bw
homer_make_ucsc_file_bigWig.U2OS_ZNF768_rep1.ZNF768.fb6fd9c14ec0650daad095d05a8fb9cd.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_6_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_6_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_7_JOB_ID: homer_make_ucsc_file.U2OS_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.U2OS_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$homer_make_tag_directory_4_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.U2OS_ZNF768_rep2.ZNF768.0679c66ea7a750422fb513d162dbc3ed.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file.U2OS_ZNF768_rep2.ZNF768.0679c66ea7a750422fb513d162dbc3ed.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.11 && \
mkdir -p tracks/U2OS_ZNF768_rep2/ZNF768 && \
touch tracks/U2OS_ZNF768_rep2/ZNF768 && \
makeUCSCfile \
  tags/U2OS_ZNF768_rep2/U2OS_ZNF768_rep2.ZNF768 \
  > tracks/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.ucsc.bedGraph && \
gzip -c tracks/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.ucsc.bedGraph \
  > tracks/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.ucsc.bedGraph.gz
homer_make_ucsc_file.U2OS_ZNF768_rep2.ZNF768.0679c66ea7a750422fb513d162dbc3ed.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_7_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 32G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_7_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_8_JOB_ID: homer_make_ucsc_file_bigWig.U2OS_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.U2OS_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$homer_make_ucsc_file_7_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.U2OS_ZNF768_rep2.ZNF768.e29386ef211693519703557925201e3c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file_bigWig.U2OS_ZNF768_rep2.ZNF768.e29386ef211693519703557925201e3c.mugqic.done' > $COMMAND
module purge && \
module load mugqic/ucsc/v346 && \
mkdir -p tracks/U2OS_ZNF768_rep2/ZNF768/bigWig && \
touch tracks/U2OS_ZNF768_rep2/ZNF768/bigWig && \
(cat tracks/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.ucsc.bedGraph | head -n 1 > tracks/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.ucsc.bedGraph.head.tmp ; ec=$?; if [ "$ec" -eq 141 ]; then exit 0; else exit "$ec"; fi) && \
cat tracks/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -v "GL\|lambda\|pUC19\|KI\|\KN\|random"  | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.ucsc.bedGraph.body.tmp && \
cat tracks/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.ucsc.bedGraph.head.tmp tracks/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.ucsc.bedGraph.body.tmp > tracks/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.ucsc.bedGraph.sorted && \
rm tracks/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.ucsc.bedGraph.head.tmp tracks/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/U2OS_ZNF768_rep2/ZNF768/bigWig/U2OS_ZNF768_rep2.ZNF768.bw
homer_make_ucsc_file_bigWig.U2OS_ZNF768_rep2.ZNF768.e29386ef211693519703557925201e3c.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_8_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_8_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_9_JOB_ID: homer_make_ucsc_file_zip
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_zip
JOB_DEPENDENCIES=$homer_make_ucsc_file_1_JOB_ID:$homer_make_ucsc_file_3_JOB_ID:$homer_make_ucsc_file_5_JOB_ID:$homer_make_ucsc_file_7_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_zip.9c17db2e15770a2905e0f8d569639e7a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file_zip.9c17db2e15770a2905e0f8d569639e7a.mugqic.done' > $COMMAND
mkdir -p report && \
zip -r report/tracks.zip tracks/*/*/*.ucsc.bedGraph.gz
homer_make_ucsc_file_zip.9c17db2e15770a2905e0f8d569639e7a.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_9_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_9_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.RAJI_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.RAJI_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$sambamba_view_filter_1_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.RAJI_ZNF768_rep1.ZNF768.a9c8aebfce737e498ca8e45e93fa9354.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak.RAJI_ZNF768_rep1.ZNF768.a9c8aebfce737e498ca8e45e93fa9354.mugqic.done' > $COMMAND
module purge && \
module load mugqic/python/3.7.3 mugqic/MACS2/2.2.7.1 && \
mkdir -p peak_call/RAJI_ZNF768_rep1/ZNF768 && \
touch peak_call/RAJI_ZNF768_rep1/ZNF768 && \
mkdir -p metrics/multiqc_inputs && \
touch metrics/multiqc_inputs && \
macs2 callpeak --format BAMPE --fix-bimodal  \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032 \
  --treatment \
  alignment/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam \
  --nolambda \
  --name peak_call/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768 \
  >& peak_call/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768.diag.macs.out && \
ln -s -f \
  ../../peak_call/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768_peaks.xls \
  metrics/multiqc_inputs/RAJI_ZNF768_rep1.ZNF768_peaks.xls
macs2_callpeak.RAJI_ZNF768_rep1.ZNF768.a9c8aebfce737e498ca8e45e93fa9354.mugqic.done
chmod 755 $COMMAND
macs2_callpeak_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 32G -c 2 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak_bigBed.RAJI_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.RAJI_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.RAJI_ZNF768_rep1.ZNF768.0f7660736cb505e0228ecb8b02f595ed.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak_bigBed.RAJI_ZNF768_rep1.ZNF768.0f7660736cb505e0228ecb8b02f595ed.mugqic.done' > $COMMAND
module purge && \
module load mugqic/ucsc/v346 && \
awk '{if ($9 > 1000) {$9 = 1000}; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)}' peak_call/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768_peaks.narrowPeak > peak_call/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/RAJI_ZNF768_rep1/ZNF768/RAJI_ZNF768_rep1.ZNF768_peaks.narrowPeak.bb
macs2_callpeak_bigBed.RAJI_ZNF768_rep1.ZNF768.0f7660736cb505e0228ecb8b02f595ed.mugqic.done
chmod 755 $COMMAND
macs2_callpeak_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.RAJI_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.RAJI_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$sambamba_view_filter_2_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.RAJI_ZNF768_rep2.ZNF768.52279289bbc424c97e5599c13fc27c63.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak.RAJI_ZNF768_rep2.ZNF768.52279289bbc424c97e5599c13fc27c63.mugqic.done' > $COMMAND
module purge && \
module load mugqic/python/3.7.3 mugqic/MACS2/2.2.7.1 && \
mkdir -p peak_call/RAJI_ZNF768_rep2/ZNF768 && \
touch peak_call/RAJI_ZNF768_rep2/ZNF768 && \
mkdir -p metrics/multiqc_inputs && \
touch metrics/multiqc_inputs && \
macs2 callpeak --format BAMPE --fix-bimodal  \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032 \
  --treatment \
  alignment/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam \
  --nolambda \
  --name peak_call/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768 \
  >& peak_call/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768.diag.macs.out && \
ln -s -f \
  ../../peak_call/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768_peaks.xls \
  metrics/multiqc_inputs/RAJI_ZNF768_rep2.ZNF768_peaks.xls
macs2_callpeak.RAJI_ZNF768_rep2.ZNF768.52279289bbc424c97e5599c13fc27c63.mugqic.done
chmod 755 $COMMAND
macs2_callpeak_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 32G -c 2 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak_bigBed.RAJI_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.RAJI_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.RAJI_ZNF768_rep2.ZNF768.4ae7b78ae0832fcad704dd0c2f8aa43a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak_bigBed.RAJI_ZNF768_rep2.ZNF768.4ae7b78ae0832fcad704dd0c2f8aa43a.mugqic.done' > $COMMAND
module purge && \
module load mugqic/ucsc/v346 && \
awk '{if ($9 > 1000) {$9 = 1000}; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)}' peak_call/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768_peaks.narrowPeak > peak_call/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/RAJI_ZNF768_rep2/ZNF768/RAJI_ZNF768_rep2.ZNF768_peaks.narrowPeak.bb
macs2_callpeak_bigBed.RAJI_ZNF768_rep2.ZNF768.4ae7b78ae0832fcad704dd0c2f8aa43a.mugqic.done
chmod 755 $COMMAND
macs2_callpeak_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.U2OS_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.U2OS_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$sambamba_view_filter_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.U2OS_ZNF768_rep1.ZNF768.dba9034384f3afebf7f296a9925318fa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak.U2OS_ZNF768_rep1.ZNF768.dba9034384f3afebf7f296a9925318fa.mugqic.done' > $COMMAND
module purge && \
module load mugqic/python/3.7.3 mugqic/MACS2/2.2.7.1 && \
mkdir -p peak_call/U2OS_ZNF768_rep1/ZNF768 && \
touch peak_call/U2OS_ZNF768_rep1/ZNF768 && \
mkdir -p metrics/multiqc_inputs && \
touch metrics/multiqc_inputs && \
macs2 callpeak --format BAMPE --fix-bimodal  \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032 \
  --treatment \
  alignment/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.sorted.dup.filtered.bam \
  --nolambda \
  --name peak_call/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768 \
  >& peak_call/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768.diag.macs.out && \
ln -s -f \
  ../../peak_call/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768_peaks.xls \
  metrics/multiqc_inputs/U2OS_ZNF768_rep1.ZNF768_peaks.xls
macs2_callpeak.U2OS_ZNF768_rep1.ZNF768.dba9034384f3afebf7f296a9925318fa.mugqic.done
chmod 755 $COMMAND
macs2_callpeak_5_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 32G -c 2 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak_bigBed.U2OS_ZNF768_rep1.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.U2OS_ZNF768_rep1.ZNF768
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.U2OS_ZNF768_rep1.ZNF768.403e84806444cda700c67d76fb10a4dc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak_bigBed.U2OS_ZNF768_rep1.ZNF768.403e84806444cda700c67d76fb10a4dc.mugqic.done' > $COMMAND
module purge && \
module load mugqic/ucsc/v346 && \
awk '{if ($9 > 1000) {$9 = 1000}; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)}' peak_call/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768_peaks.narrowPeak > peak_call/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/U2OS_ZNF768_rep1/ZNF768/U2OS_ZNF768_rep1.ZNF768_peaks.narrowPeak.bb
macs2_callpeak_bigBed.U2OS_ZNF768_rep1.ZNF768.403e84806444cda700c67d76fb10a4dc.mugqic.done
chmod 755 $COMMAND
macs2_callpeak_6_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$macs2_callpeak_6_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.U2OS_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.U2OS_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$sambamba_view_filter_4_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.U2OS_ZNF768_rep2.ZNF768.bdb047ac3012f642a662092e33ae7af8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak.U2OS_ZNF768_rep2.ZNF768.bdb047ac3012f642a662092e33ae7af8.mugqic.done' > $COMMAND
module purge && \
module load mugqic/python/3.7.3 mugqic/MACS2/2.2.7.1 && \
mkdir -p peak_call/U2OS_ZNF768_rep2/ZNF768 && \
touch peak_call/U2OS_ZNF768_rep2/ZNF768 && \
mkdir -p metrics/multiqc_inputs && \
touch metrics/multiqc_inputs && \
macs2 callpeak --format BAMPE --fix-bimodal  \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032 \
  --treatment \
  alignment/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.sorted.dup.filtered.bam \
  --nolambda \
  --name peak_call/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768 \
  >& peak_call/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768.diag.macs.out && \
ln -s -f \
  ../../peak_call/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768_peaks.xls \
  metrics/multiqc_inputs/U2OS_ZNF768_rep2.ZNF768_peaks.xls
macs2_callpeak.U2OS_ZNF768_rep2.ZNF768.bdb047ac3012f642a662092e33ae7af8.mugqic.done
chmod 755 $COMMAND
macs2_callpeak_7_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 32G -c 2 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$macs2_callpeak_7_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak_bigBed.U2OS_ZNF768_rep2.ZNF768
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.U2OS_ZNF768_rep2.ZNF768
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.U2OS_ZNF768_rep2.ZNF768.e30070f92a7f347ba37437954b63431b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak_bigBed.U2OS_ZNF768_rep2.ZNF768.e30070f92a7f347ba37437954b63431b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/ucsc/v346 && \
awk '{if ($9 > 1000) {$9 = 1000}; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)}' peak_call/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768_peaks.narrowPeak > peak_call/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/root/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/U2OS_ZNF768_rep2/ZNF768/U2OS_ZNF768_rep2.ZNF768_peaks.narrowPeak.bb
macs2_callpeak_bigBed.U2OS_ZNF768_rep2.ZNF768.e30070f92a7f347ba37437954b63431b.mugqic.done
chmod 755 $COMMAND
macs2_callpeak_8_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
GenPipes_STATE=\$PIPESTATUS
echo GenPipesExitStatus:\$GenPipes_STATE



if [ \$GenPipes_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (GenPipes)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$GenPipes_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 3900M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$macs2_callpeak_8_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
LOG_MD5=$(echo $USER-'10.80.49.1-ChipSeq-RAJI_ZNF768_rep1.ChIP_seq_RAJI_ZNF768_rep1,RAJI_ZNF768_rep2.ChIP_seq_RAJI_ZNF768_rep2,U2OS_ZNF768_rep1.ChIP_seq_U2OS_ZNF768_rep1,U2OS_ZNF768_rep2.ChIP_seq_U2OS_ZNF768_rep2' | md5sum | awk '{ print $1 }')
if test -t 1; then ncolors=$(tput colors); if test -n "$ncolors" && test $ncolors -ge 8; then bold="$(tput bold)"; normal="$(tput sgr0)"; yellow="$(tput setaf 3)"; fi; fi
wget --quiet 'http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=narval1.narval.calcul.quebec&ip=10.80.49.1&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,mapping_bwa_mem_sambamba,sambamba_merge_bam_files,sambamba_mark_duplicates,sambamba_view_filter,bedtools_blacklist_filter,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak&samples=4&md5=$LOG_MD5' -O /dev/null || echo "${bold}${yellow}Warning:${normal}${yellow} Genpipes ran successfully but was not send telemetry to mugqic.hpc.mcgill.ca. This error will not affect genpipes jobs you have submitted.${normal}"
