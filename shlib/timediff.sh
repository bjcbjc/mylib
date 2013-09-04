#!/bin/bash


ext1=$1
ext2=$2

for f in $(ls *.$ext1); do
    fn="${f%.*}"
    if [ -e $fn.$ext2 ]; then
	hr=$(echo "scale=1; ($(date -r $fn.$ext2 +%s) - $(date -r $f +%s))/3600 " | bc)
	echo -e "$hr hrs\t$fn"
    fi
done


