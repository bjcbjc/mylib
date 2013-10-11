#!/bin/bash
#$ -cwd
#$ -j y

bam=$1
if [ $# -lt 2 ]; then
    logfile=$bam".split.log"
else
    logfile=$2
fi

chrmformat=`/data/NYGC/Software/samtools/samtools-0.1.19/samtools view $bam | head -1 | cut -f3`
if [[ "$chrmformat" == chr* ]]; then
    nfile=`wc -l /nethome/bjchen/DATA/GenomicInfo/chrmsplit.forbam_chrtag | cut -f1 -d" "`
else
    nfile=`wc -l /nethome/bjchen/DATA/GenomicInfo/chrmsplit.forbam | cut -f1 -d" "`
fi

qsub -t 1-$nfile -o $logfile /nethome/bjchen/BJLib/shlib/splitbam.subscript.sh $bam


