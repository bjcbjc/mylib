

infile=$1
outfile=${infile/.bam/.AS.bam}
outfile=${outfile#/*bam}
inpath=${infile%/*bam}
outpath=/scratch/BJC_TMP/$RANDOM

mkdir $outpath

/data/NYGC/Software/samtools/samtools-0.1.19/samtools view -h $infile |  gawk '{OFS="\t"; if($0 ~"^@") {print $0} else{ for(i=1;i<=NF;i++){if($i ~/AS:i/) { $5=$i; sub(/AS:i:/,"",$5); print $0}}}}' | /data/NYGC/Software/samtools/samtools-0.1.19/samtools view -Sbh - > $outpath$outfile

mv $outpath$outfile $inpath$outfile
rmdir $outpath

