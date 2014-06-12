#extract single sample, filter 0/0 according to the sample
#gawk -v sample="" -f vcf.getsample.awk vcf_file

BEGIN { count = 0 }
$1 ~ /^##/ { \
    header[++count] = $0
}
$1 ~ /^#CHROM/ { \
    sampleidx=0
    for (i=1; i<=NF; i++) {
	if ($i == sample) { 
	    sampleidx = i
	    break
	}
    }
    if (sampleidx == 0) {
	print "Cannot find sample "sample
	exit 1
    }else {
	for (i=1; i<=length(header); i++) { print header[i] }
	delete header
	for (i=1; i<=9; i++) { printf("%s\t",$i) }
	printf("%s\n",$sampleidx)
    }
}
$1 !~ /^#/ { \
    if (index($sampleidx, "0/0") == 0 && index($sampleidx, "./.") == 0 ) {
	for (i=1; i<=9; i++) {
	    printf("%s\t",$i)
	}
	printf("%s\n",$sampleidx)
    }
}

