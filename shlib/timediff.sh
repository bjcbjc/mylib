#!/bin/bash


ext1=$1
ext2=$2
if [ $# -lt 3 ]; then
    dir=""
else
    dir=$3
fi


for f in $(ls $dir*.$ext1); do
    fn="${f%.*}"
    if [ -e $fn.$ext2 ]; then
	hr=$(echo "scale=1; ($(date -r $fn.$ext2 +%s) - $(date -r $f +%s))/3600 " | bc)
	echo -e "$hr hrs\t$fn"
    fi
done


