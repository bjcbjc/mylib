
$1 !~ /^#/ { \
    OFS="\t"
    indel=0
    for (col=4; col<=5; col++) {
	if ($col ~ /,/) {
	    split($col, tmp, ",")
	    for (i=1; i<=length(tmp); i++) {
		if (length(tmp[i]) > 1) {
		    indel=1
		    break
		}
	    }
	} else {
	    if (length($col) > 1) {
		indel=1
	    }
	}
	if (indel==1)
	    break
    }
    if (indel==1)
	print $1,$2,$4,$5
}

