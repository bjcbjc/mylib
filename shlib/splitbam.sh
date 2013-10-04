#!/bin/bash
#$ -cwd
#$ -j y

bam=$1
if [ $# -lt 2 ]; then
    logfile=$bam".split.log"
else
    logfile=$2
fi

nfile=`wc -l /nethome/bjchen/DATA/GenomicInfo/chrmsplit.forbam | cut -f1 -d" "`

qsub -t 1-$nfile -o $logfile /nethome/bjchen/BJLib/shlib/splitbam.subscript.sh $bam


