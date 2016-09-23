
import sys
import subprocess
import os.path

old = sys.argv[1]
new = sys.argv[2]

tb = dict()
for line in open(new):
    line = line.strip().split()
    if len(line) != 3: continue
    jobname = line[0]
    tb[jobname] = [line[1], False]

count = 0
data = [line.strip().split('\t') for line in open(old)]
## copy file in case of errors
if os.path.exists(old + '.copy'):
    print '%s.copy already exists. Cannot proceed.'%(old)
    exit()

error = subprocess.call(['cp', old, old + '.copy'])
if error != 0:
    print 'Cannot make a copy of the original job id: %s'%error
    exit()

restore = False
try:    
    with open(old, 'w') as f:
        for line in data:
            job = line[0]
            if job in tb:
                tb[job][1] = True
                line[1] = tb[job][0]
                count += 1
            f.write('\t'.join(line) + '\n')
except:
    restore = True

notReplaced = sum([not val[1] for val in tb.values()])
if count != len(tb) or notReplaced != 0:
    restore = True
    print 'error: replaced %d, len(tb)=%d, notReplaced:%d, revert the original (no changes made)'%(count, len(tb), notReplaced)
else:
    print 'replaced %d'%count

if restore:
    subprocess.call(['mv', old + '.copy', old])
else:
    subprocess.call(['rm', '-f', old + '.copy'])
