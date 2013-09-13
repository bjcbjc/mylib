#!/bin/bash

# simplejob.sh   "command"  "file patterns used in ls"
#
# Creates a temporary script for SGE array jobs; one job-ID is associated
# with all jobs where each execute "command" with one of the files. 
#

cmd=$1
filepattern=$2
runmode=$3

if [ $# -lt 4 ]; then
    pipeoutext=""
else
    pipeoutext=$4
    if [[ "$pipeoutext" != .* ]]; then
	pipeoutext="."$pipeoutext
    fi
fi

tmpstr=$RANDOM

ls $filepattern > simpjob.$tmpstr.input
nfile=`wc -l simpjob.$tmpstr.input | cut -f1 -d" "`
echo "list in simpjob.$tmpstr.input"

tmpscript=simpjob.$tmpstr.job
echo "#!/bin/bash" > $tmpscript
echo "#$ -cwd" >> $tmpscript
echo "#$ -j y" >> $tmpscript
echo "#$ -o simpjob.$tmpstr.out" >> $tmpscript
echo "infile=\`awk \"NR==\$SGE_TASK_ID\" simpjob.$tmpstr.input\`" >> $tmpscript

if [ -z "$pipeoutext" ]; then
    echo "$cmd \$infile " >> $tmpscript
else
    echo "$cmd \$infile > \$infile$pipeoutext" >> $tmpscript
fi

echo "echo task \$SGE_TASK_ID done" >> $tmpscript
if [ $runmode = "run" ]; then    
    qsub -t 1-$nfile $tmpscript
fi




