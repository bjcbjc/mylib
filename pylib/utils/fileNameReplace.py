#!/usr/bin/python

from os import popen, system
from sys import argv, exit
import os.path
import re

if len(argv) < 3:
    print 'fileNameReplace.py <files> <string to replace> <new string>'
    exit()

f = popen('ls %s'%argv[1])
fns = re.split('\t|\n', f.read())
f.close()

for a in fns:
    a = a.replace('(','\(').replace(')','\)').replace(',','\,').replace(' ','\ ')
    #a = re.sub('(\W)', '\\\0',a)
    newname = a.replace(argv[2], argv[3])
    if os.path.exists(newname):
        print '%s exists, skip %s -> %s'%(newname, a, newname)
    else:
        system('mv %s %s'%(a, newname))

