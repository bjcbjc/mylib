
vcf=$1
grep -v ^# $vcf | gawk '{if(length($4)==1 && length($5)==1) print $4"->"$5 }' | sort | uniq -c
