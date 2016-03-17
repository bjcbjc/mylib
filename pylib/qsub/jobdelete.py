#!/usr/bin/python

from os import popen, system
from sys import argv

keyword = argv[1]
#f = popen('qstat -u bjchen | grep %s'%keyword)
qstat2='/data/NYGC/Software/java/jre1.7.0_25/bin/java -jar /data/NYGC/software/SGETools/sgetools.jar -a qstat'
f = popen('%s -u bjchen | grep %s'%(qstat2, keyword))
t = f.read()
f.close()

t = t.strip().split('\n')
jid = map(lambda(l):l.split()[0],t)
print 'going to delete %d jobs related %s'%(len(jid),keyword)
for id in jid:
    system('qdel %s'%id)
