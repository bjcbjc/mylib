BEGIN {
    OFS="\t"
}
$0 ~ /^#/ {
    print $0
}
$0 !~ /^#/ {					
    #MT
    if ($1 == "MT") {
	$1 = "chrM"
    } else if ($1 ~ /[0-9]+/) {
	$1 = "chr"$1
    }
    print $0
}

