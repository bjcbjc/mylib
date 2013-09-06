#!/bin/bash

cmd=$1
filepattern=$2

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
echo "$cmd \$infile" >> $tmpscript
echo "echo task \$SGE_TASK_ID done" >> $tmpscript
qsub -t 1-$nfile $tmpscript




