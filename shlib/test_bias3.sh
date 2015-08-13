#!/bin/bash

if [ $# -ne 4 ]; then
    echo "Usage: bash test_bias3.sh exon_count_file output_path output_prefix gene_body_cov"
    exit 1
fi


set -e

Rcmd=/data/NYGC/Software/R/R-3.1.0/bin/Rscript
script=/nethome/bjchen/Projects/bench/genebodycov/rscript/test3PrimeBias.r
exonfn=${1}
outputDir=${2}
output=${3}
gbfn=${4}
log=${outputDir}/${output}.Exon3BiasTest.log

## set up 
if [ ! -d "${outputDir}" ]; then 
mkdir -p "${outputDir}"
fi


cmd="qsub -l h_vmem=10G -b y -j y -o ${log} -N ${output} ${Rcmd} ${script} input=${exonfn} outputPath=${outputDir} outputPrefix=${output} gbfn=${gbfn}"
echo ${cmd}
${cmd}

