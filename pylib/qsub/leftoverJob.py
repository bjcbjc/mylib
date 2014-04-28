
from os import popen
from sys import argv

#usage" leftoverJob.py dir_with_jobs
#
# for all the jobs in the directory, see which jobs are still in queue and return those jobs that are not in queue

def getLeftoverJob( jobdir ):

    cmd = '/data/NYGC/Software/java/jre1.7.0_25/bin/java -jar /data/NYGC/software/SGETools/sgetools.jar -a qstat'

    f = popen(cmd)
    qjob = [ line.split()[1] for line in f.readlines()[2:-1] ]
    f.close()

    try:
        f = popen('ls %s/*.job'%( jobdir ) )
        job = [ line.split('/')[-1] for line in f.read().split() ]
        f.close()
    except: 
        job = []

    if len(job) == 0:
        print 'no job in %s'%jobdir
    leftover = list( set(job).difference(qjob) )
    return leftover, job, qjob


if __name__ == '__main__':
    if len(argv) == 1:
        jobdir = './'
    else:
        jobdir = argv[1]
    leftover, job, qjob = getLeftoverJob( jobdir )    
    if len(leftover) > 0:
        for line in leftover:
            print line
    else:
        print '%d jobs in the directory; %d in queue'%(len(job), len(qjob))

