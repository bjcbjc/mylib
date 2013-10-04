
bam=$1
if [ $# -lt 2 ]; then
    n=10
else
    n=$2
fi

flist=`eval echo $bam.{2..$n}`

/data/NYGC/Software/samtools/samtools-0.1.19/samtools view -H $bam | tee $flist >/dev/null

/data/NYGC/Software/samtools/samtools-0.1.19/samtools view -h $bam | gawk -v nfile="$n" -v bamroot="$bam" '{if(NR < 500){ print $0 > bamroot".1"} else{part=(NR-1)%nfile+1; print $0 >> bamroot"."part} }'
#samtools view -b $bam | awk -v nfile="$n" -v bamroot="$bam" '{if (NR < 500) print $0 >> bamroot".1"; else {part=(NR-1)%nfile+1; print $0 >> bamroot"."part} }'
