#!/bin/bash

for i in {1..31}; do
    node=`printf node%03d $i`
    echo $node
    ssh $node bash /nethome/bjchen/BJLib/shlib/rmEmptyDirRecur.sh /scratch/BJC_TMP/
done
