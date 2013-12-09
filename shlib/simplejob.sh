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
    jobname=$RANDOM
else
    jobname=$4
fi

if [ $# -lt 5 ]; then
    pipeoutext=""
else
    pipeoutext=$5
    if [[ "$pipeoutext" != .* ]]; then
	pipeoutext="."$pipeoutext
    fi
fi

ls $filepattern > simpjob.$jobname.input
nfile=`wc -l simpjob.$jobname.input | cut -f1 -d" "`
echo "list in simpjob.$jobname.input"

tmpscript=simpjob.$jobname.job
echo "#!/bin/bash" > $tmpscript
echo "#$ -cwd" >> $tmpscript
echo "#$ -j y" >> $tmpscript
echo "#$ -o simpjob.$jobname.out" >> $tmpscript
echo "infile=\`awk \"NR==\$SGE_TASK_ID\" simpjob.$jobname.input\`" >> $tmpscript

if [ -z "$pipeoutext" ]; then
    if [[ "$cmd" =~ "\$infile" ]]; then
	echo "$cmd " >> $tmpscript
    else
	echo "$cmd \$infile " >> $tmpscript
    fi
else
    if [[ "$cmd" =~ "\$infile" ]]; then
	echo "$cmd > \$infile$pipeoutext" >> $tmpscript
    else
	echo "$cmd \$infile > \$infile$pipeoutext" >> $tmpscript
    fi
fi

echo "echo task \$SGE_TASK_ID done" >> $tmpscript
if [ $runmode = "run" ]; then    
    qsub -t 1-$nfile $tmpscript
fi




