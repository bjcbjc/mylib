#!/bin/bash
#$ -cwd
#$ -j y


infile=`awk "NR==$SGE_TASK_ID" /nethome/bjchen/DATA/GenomicInfo/chrmsplit.forbam`

echo "/data/NYGC/Software/samtools/samtools-0.1.19/samtools view -bh $1 $infile > "$1"."$SGE_TASK_ID".bam"
/data/NYGC/Software/samtools/samtools-0.1.19/samtools view -bh $1 $infile > $1"."$SGE_TASK_ID".bam"

echo task $SGE_TASK_ID done
