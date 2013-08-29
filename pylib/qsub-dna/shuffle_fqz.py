


from subprocess import Popen, PIPE
from sys import argv

file1 = argv[1]
file2 = argv[2]
out = argv[3]

f1 = Popen('zcat %s'%file1, shell=True, stdout=PIPE).stdout
f2 = Popen('zcat %s'%file2, shell=True, stdout=PIPE).stdout
o = Popen('gzip - > %s'%out, shell=True, stdin=PIPE).stdin

l = f1.readline()
while l:
    o.write('%s'%l)
    for i in range(3):
        o.write(f1.readline())
    for i in range(4):
        o.write(f2.readline())
    l = f1.readline() #next read?

f1.close()
f2.close()
o.close()
    
