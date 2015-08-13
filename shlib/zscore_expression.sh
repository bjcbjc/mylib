#!/bin/bash

if [ $# -lt 3 ]; then
    echo "Usage: bash zscore_expression.sh  cohort_fn count_fn cohortName [valid_cohort_sample_list]"
    exit 1
fi


set -e

matlab=/usr/local/MATLAB/R2014a/bin/matlab
cohort=${1}
rnaCount=${2}
cohortName=${3}
validSampleList=${4}

## set up 
output=${rnaCount///qc//expression}
output=${output///featureCounts//expression}
output=${output/_counts.txt/}
outputDir=$(dirname ${output})
log=${output}.${cohortName}.log

if [ ! -d "${outputDir}" ]; then 
mkdir -p "${outputDir}"
fi


exec &>${log}
echo hostname: `hostname`
hostname=`hostname`


echo "===== normalize expression : "$(date)

matlabcall="${matlab} -nodisplay -r \"try; addpath(genpath('/nethome/bjchen/BJLib/Matlabox')); zscore_expression('${cohort}', '${rnaCount}', '${output}', '${cohortName}', '${validSampleList}'); fprintf('Done normalization\n'); quit; catch err; fprintf('error\n'); rethrow(err); exit(1); end\""

echo ${matlabcall}

if [ $hostname = "matlab.nygenome.org" ]; then
	${matlabcall}
else
	ssh matlab.nygenome.org "exec &>>${log}; ${matlabcall}"
fi


