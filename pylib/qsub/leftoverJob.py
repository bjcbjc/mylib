
from os import popen
from sys import argv

#usage" leftoverJob.py dir_with_jobs
#
# for all the jobs in the directory, see which jobs are still in queue and return those jobs that are not in queue

cmd = '/data/NYGC/Software/java/jre1.7.0_25/bin/java -jar /data/NYGC/software/SGETools/sgetools.jar -a qstat'

f = popen(cmd)
qjob = [ line.split()[1] for line in f.readlines()[2:-1] ]
f.close()

f = popen('ls %s/*.job'%(argv[1]))
job = [ line.split('/')[-1] for line in f.read().split() ]
f.close()

leftover = list( set(job).difference(qjob) )

if len(leftover) > 0:
    for line in leftover:
        print line
else:
    print '%d jobs in the directory; %d in queue'%(len(job), len(qjob))



