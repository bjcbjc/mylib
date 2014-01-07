
import argparse
import subprocess
import traceback
import time
from collections import Counter

samflaglist = [ (0x1, 'read paired'), (0x2, 'read mapped in proper pair'), (0x4, 'read unmapped'), (0x8, 'mate unmapped'), \
              (0x10, 'read is reverse-complemented'), (0x20, 'mate is reverse-complemented'), (0x40, 'first in pair'), \
              (0x80, 'second in pair'), (0x100, 'secondary alignment'), (0x200, 'failed platform QC'), \
              (0x400, 'PCR/optical duplicate'), (0x800, 'supplementary aglignment') ]
v, k = zip(*samflaglist)
samflagtb = dict(zip(k,v))
summarizeFlag = {'read mapped in proper pair': 0x2, 'read unmapped': 0x4, 'secondary alignment': 0x100}


# helper functions
def cumsum(l):
    s = l[:] #make a copy
    for i in xrange(1, len(s)):
        s[i] = s[i]+s[i-1]
    return s

def closefiles(files):
    for f in files:
        if f != None: f.close()
    return

# set up parameters for the script
argp = argparse.ArgumentParser(prog='bamFlagStat.py')
argp.add_argument('-bam', type=str, default = '', metavar='bam_file', required=False, help='input file')
argp.add_argument('-flagcount', type=str, default = '', metavar='flag_count_file', required=False, help='input file')
argp.add_argument('-o', type=str, default='bam.flag.stat', metavar='output_file', help='file name for the output, [<bam>.flag.stat]')
argp.add_argument('-samtools', type=str, default='/data/NYGC/Software/samtools/samtools-0.1.19/samtools', help='samtools path, [/data/NYGC/Software/samtools/samtools-0.1.19]')

args = argp.parse_known_args()
args = args[0]

hasbam = len(args.bam) > 0
hasflagcount = len(args.flagcount) > 0

if hasflagcount:
    if hasbam:
        print 'ignore bam file since flag-count file is provided'
    flagcount = [ map(int, line.strip().split()) for line in open(args.flagcount).readlines() ] #count, flag
else:
    if hasbam:
        print 'do flag-count first'
        cmd = args.samtools + 'view ' + args.bam + ' | cut -f2 | sort | uniq -c '
        try:
            samout = None
            sampipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            samout = sampipe.stdout
            flagcount = samout.read()
            samout.close()
        except:
            closefiles( [samout] )
            traceback.print_exc()
        f = open('%s.flagcount'%args.o, 'w')
        f.write('%s'%flagcount)
        f.close()
        flagcount = [ map(int, line.strip().split()) for line in flagcount.strip().split('\n') ]  #count, flag
    else:
        print 'need to provide either bam or flag-count file; flag-count file has two columns: the first contains the count of alignments and the second the flag (integer) in the bam'
        exit()

translation = []
summary = Counter()
for count, flag in flagcount:
    translation.append( [count, ', '.join( [ desc for b, desc in samflaglist if (flag & b) > 0 ] ) ] )
    for sumkey in summarizeFlag.keys():
        if (flag & summarizeFlag[sumkey]) > 0:
            summary[sumkey] = summary[sumkey] + count

translation = sorted( translation, key=lambda line: line[1] )

f = open(args.o, 'w')
f.write('===== Summary =====\n')
for sumkey in summarizeFlag.keys():
    f.write('%10d\t\t%s\n'%(summary[sumkey], sumkey))
f.write('===== Counts of all alignments =====\n')
for count, desc in translation:
    f.write('%10d\t\t%s\n'%(count, desc))
f.close()


    
