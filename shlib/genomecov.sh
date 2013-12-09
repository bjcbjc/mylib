
bedtools_version="bedtools-2.17.0"
genome="hg19"
infile=$1
outfile=${infile/.bam/.genome.cov}
outfile=${outfile#/*genome.cov}
inpath=${infile%/*bam}
outpath=/scratch/BJC_TMP/$RANDOM


mkdir $outpath

/data/NYGC/Software/bedtools/$bedtools_version/bin/bedtools genomecov -bg -split -ibam $infile -g /data/NYGC/Software/bedtools/$bedtools_version/genomes/human.$genome.genome > $outpath$outfile

mv $outpath$outfile $inpath$outfile
rmdir $outpath

