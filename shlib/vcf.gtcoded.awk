## parse vcf to chrm, pos, type, filter, gt
BEGIN { 
    OFS = "\t"
    if (length(prefix) == 0) {
	print "need prefix"
    }
    gtTable["./."] = -1
    gtTable["0/0"] = 0
    gtTable["0/1"] = 1
    gtTable["1/1"] = 2
    nextGtCode = 3
    maskFn = sprintf("%s.GT.mask.txt", prefix)
    codeFn = sprintf("%s.GT.codes.txt", prefix)
    if (length(info) > 0) {
	split(info, infoFields, ",")
    }
}

$1 == "#CHROM" {
    printf("%s\t%s\t%s\t%s\t%s\tTYPE", $1, $2, $4, $5, $7) > maskFn
    for (i = 1; i <= length(infoFields); i++) {
	printf("\t%s", infoFields[i]) > maskFn
    }
    for (i = 10; i <= NF; i++) {
	printf("\t%s", $i) > maskFn
    }
    printf("\n") > maskFn
}

$1 !~ /^#/ {
    #1:chr, 2:pos, 4,5:ref/alt, 7:filter, 9:format 10-: sample
    delete genotype
    delete infoData

    ## find info
    for (i=1; i <= length(infoFields); i++) {
       found = match($8, infoFields[i]"=([a-zA-Z0-9\\-\\+\\.\\,\\_]+)", tmpInfoData)
       if (found != 0) {
	   infoData[i] = tmpInfoData[1]
       } else { 
	   infoData[i] = "NA"
       }
   }
    
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
	    gtStr = unphase(data[gtIdx])
	    if (!(gtStr in gtTable)) {
		gtTable[gtStr] = nextGtCode
		nextGtCode += 1
	    }
	    genotype[colIdx - 9] = gtTable[gtStr]
	}
    }

    vtype = var_type($4, $5)

    ## output $1, $2, $4, $5, $7, vtype, genotype
    ##printf("%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $7, $4, $5, vtype)
    printf("%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $4, $5, $7, vtype) > maskFn

    ## info
    for (i = 1; i <= length(infoData); i++)
	printf("\t%s", infoData[i]) > maskFn
    
    ## GT
    for (i = 1; i <= length(genotype); i++) 
	printf("\t%s", genotype[i]) > maskFn
    printf("\n") > maskFn
}

END {
    for (gtStr in gtTable) 
	printf("%s\t%s\n", gtStr, gtTable[gtStr]) > codeFn
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

function unphase(gtStr) {
    if (index(gtStr, "|") == 0) {
	split(gtStr, data, "/")
	if (data[1] < data[2])
	    return data[1]"/"data[2]
	else
	    return data[2]"/"data[1]
    } else {
	split(gtStr, data, "|")
	if (data[1] < data[2])
	    return data[1]"/"data[2]
	else
	    return data[2]"/"data[1]
    }
}
