## parse vcf to chrm, pos, type, filter, gt
BEGIN { 
    OFS = "\t"
#    if (length(prefix) == 0) {
#	print "need prefix"
#    }
}

$1 == "#CHROM" {
    printf("%s\t%s\t%s\t%s\t%s\tTYPE", $1, $2, $4, $5, $7)
    for (i = 10; i <= NF; i++) {
	printf("\t%s", $i)
    }
    printf("\n")
}

$1 !~ /^#/ {
    #1:chr, 2:pos, 4,5:ref/alt, 7:filter, 9:format 10-: sample
    delete genotype
    
    ## find GT
    split($9, format, ":")
    gtIdx = 0
    for (i = 1; i <= length(format); i++ ) {
	if (format[i] == "GT") {
	    gtIdx = i
	}
    }
    
    for (colIdx = 10; colIdx <= NF; colIdx++) {
	genotype[colIdx - 9] = "NA"  #initialize just in case we can't find GT
	if ( gtIdx > 0) {
	    split($colIdx, data, ":")
	    genotype[colIdx - 9] = data[gtIdx]
	}
    }

    vtype = var_type($4, $5)

    ## output $1, $2, $4, $5, $7, vtype, genotype
    ##printf("%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $7, $4, $5, vtype)
    printf("%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $4, $5, $7, vtype)
    for (i = 1; i <= length(genotype); i++) 
	printf("\t%s", genotype[i])
    printf("\n")

}

function var_type(ref, alt) {
    refLen = length(ref)
    altLen = length(alt)
    multiRef = index(ref, ",")
    multiAlt = index(alt, ",")
    if (multiRef > 0 || multiAlt > 0) {
	vtype = "MULTI"
    } else if (refLen != altLen) {
	vtype = "INDEL"
    } else if (refLen == 1 && altLen == 1) {
	vtype = "SNP"
    } else {
	vtype = "OTHER"
    }
    return vtype
}
