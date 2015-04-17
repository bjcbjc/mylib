
BEGIN { OFS="\t" }


FNR > 0 {
    if (length(_[FNR]))
	_[FNR]=(_[FNR] OFS $2)
    else
	_[FNR]=$0
}

END {
    for(x=1; x <= length(_); x++)
	if (length(_[x]))
	    print _[x]
}
