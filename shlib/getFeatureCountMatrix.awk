
BEGIN { OFS="\t" }

FNR == 2 {
    match($7, /Sample_([0-9A-Za-z\-\_]+)/, sample)
    $7 = sample[1]
}
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
