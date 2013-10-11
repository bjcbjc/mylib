#!/bin/bash
#$ -cwd
#$ -j y

bam=$1
chrmformat=`/data/NYGC/Software/samtools/samtools-0.1.19/samtools view $bam | head -1 | cut -f3`
if [[ "$chrmformat" == chr* ]]; then
    infile=`awk "NR==$SGE_TASK_ID" /nethome/bjchen/DATA/GenomicInfo/chrmsplit.forbam_chrtag`
else
    infile=`awk "NR==$SGE_TASK_ID" /nethome/bjchen/DATA/GenomicInfo/chrmsplit.forbam`
fi



echo "/data/NYGC/Software/samtools/samtools-0.1.19/samtools view -bh $1 $infile > "$1"."$SGE_TASK_ID".bam"
/data/NYGC/Software/samtools/samtools-0.1.19/samtools view -bh $1 $infile > $1"."$SGE_TASK_ID".bam"

echo task $SGE_TASK_ID done
