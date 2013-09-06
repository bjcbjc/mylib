#!/bin/bash


file=$1
field=$2
repstr=$3

if [ $# -lt 4 ]; then
    sampara=""
else
    sampara=$4
fi

samtools view $sampara -h $file | awk '/^@/ {print}
/^[^@]/ {$'$field'="'$repstr'"; OFS="\t"; print}' - | samtools view -Sb -
