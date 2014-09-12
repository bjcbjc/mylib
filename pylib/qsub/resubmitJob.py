
#some jobs failed because of permission issues
#this script is to resubmit these jobs without calling runjob.py again
#
#

from os import system
from sys import argv
from leftoverJob import getLeftoverJob

jobdir = argv[1]
jobfn, jobInDir, jobInQ = getLeftoverJob( jobdir )

njob = len(jobfn)

while 1:
    if njob <= 10:
        proceed = raw_input('Are you sure to resubmit %d jobs %s (yes/no): '%(njob, jobfn) )
    else:
        proceed = raw_input('Are you sure to resubmit %d jobs (yes/no): '%(njob) )
    if proceed in ['yes', 'no']:
        break


if proceed == 'yes':
    for fn in jobfn:
        #clean the sample directory and make sure the log directory exists
        t = open(fn).readlines()
        logfn = [ line for line in t if '#$ -o ' in line][0].replace('#$ -o ', '')
        sampledir = logfn.split('/logs')[0] + '/'
        system('rm -f %s'%logfn)
#        system('rm -rf %s*'%sampledir)
#        system('mkdir -p %slogs'%sampledir)
        system('qsub %s'%fn)
    
