
BEGIN { OFS="\t" }

FNR > 1 {
    if (length(_[FNR]))
	_[FNR]=(_[FNR] OFS $7)
    else
	_[FNR]=$0
}

END {
    for(x=1; x <= length(_); x++)
	if (length(_[x]))
	    print _[x]
}
