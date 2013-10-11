

import argparse
import subprocess
import traceback
import re
import time
import gzip
from collections import Counter
from sys import exit

#output:
#chr, pos, total_read, [num_read_with_var_at_1st_read_base, ...]*4 (R1F,R1R,R2F,R2R)



# helper functions
def cumsum(l):
    s = l[:] #make a copy
    for i in xrange(1, len(s)):
        s[i] = s[i]+s[i-1]
    return s

def outputline(chrm, pos, numread, counter, readlength):
    keys = ['first-fwd', 'first-rev', 'second-fwd', 'second-rev']
    line = '%s\t%s\t%s'%(chrm, pos, numread)
    for k in keys:
        for i in xrange(1, readlength+1):
            line = line + '\t%d'%(counter[k][i])
    line = line + '\n'
    return line

def closefiles(files):
    for f in files:
        if f != None: f.close()
    return

# set up parameters for the script
argp = argparse.ArgumentParser(prog='varPositionInRead.py')
argp.add_argument('-bam', type=str, metavar='bam_file', required=True, help='bam file')
argp.add_argument('-o', type=str, default='test.var.pos.gz', metavar='output_file', help='file name for the output')
argp.add_argument('-vcf', type=str, metavar='vcf_file', required=True, help='vcf file')
argp.add_argument('-q', type=int, metavar='min_mapq', help='min MAPQ', default=1)
argp.add_argument('-r', type=int, metavar='read_length', default=0, help='max read length; if not provided, obtained from the first read')
argp.add_argument('-samtools', type=str, metavar='path_to_samtools',default='/data/NYGC/Software/samtools/samtools-0.1.19/samtools', help='samtools path, [/data/NYGC/Software/samtools/samtools-0.1.19]')

args = argp.parse_known_args()
passon = ''
if isinstance(args, tuple):
    passon = ' '.join(args[1])
    args = args[0]
print args

if args.o[-3:] != '.gz' and args.o[-5:] != '.gzip':
    args.o = args.o + '.gz'

# Command string; "samtools view -F" is used to filter reads and retain only flag, pos, CIGAR, and seq. 
cmd = args.samtools + ' view -f 0x2 -F %d -q %d '%(0xf04, args.q) + args.bam + ' %s:%s-%s ' + ''' | gawk '{OFS="\t"; print $2,$4,$6,$10}' '''
#cmd = args.samtools + ' view ' + args.bam + ' %s:%s-%s ' + ''' | gawk '{OFS="\t"; print $2,$4,$6,$10}' '''

#print 'Run command: ' + cmd

# set up pipes
samout, samerr, out = None, None, None
vcf = open(args.vcf)

# container for statistics to be collected
counter = {'first-fwd':Counter(), 'first-rev':Counter(), 'second-fwd':Counter(), 'second-rev':Counter()}
# regular expression parser for CIGAR and MD
cigarparser = re.compile('(\d+)([MIDXSN\=])')

#in CIGAR, sum of M, I, X, = and S equals to read length; N (gap such as intron) and D (deletion) are not in reads
#to calculate read-base's mapped genomic location, we need to include N and D but not I, since insertions are not in reference
#inserted read-base can take n+0.5 where n is the leftmost (5') genomic location where insertion occurs

start = time.time()
try:
    # decide read length first
    if args.r == 0:
        sampipe = subprocess.Popen(args.samtools + ' view ' + args.bam , shell=True, stdout=subprocess.PIPE)
        samout = sampipe.stdout #, sampipe.stderr
        readlength = len( samout.readline().split()[9] )
        samout.close()
    else:
        readlength = args.r

    out = gzip.open(args.o, 'wb')

    nvariant = 0
    outputbuffer = ''
    for varline in vcf:
        if varline[0] == '#': continue
        varline = varline.split()
        varchrm, varpos, alt = varline[0], varline[1], varline[4]

        #implement later
        if ',' in alt: alt = alt.split(',')[0]
        if len(alt) > 1: continue

        nvariant = nvariant + 1
        #reset counters
        for c in counter.values(): c.clear()

        sampipe = subprocess.Popen(cmd%(varchrm, varpos, varpos), shell=True, stdout=subprocess.PIPE)
        samout = sampipe.stdout #, sampipe.stderr

        varpos = int(varpos)
        line = samout.readline()
        numread = 0
        while line:
            numread  = numread + 1
            flag, pos, cigar, seq = line.split()
            flag = int(flag)
            pos = int(pos)

            if cigar == '*': #not available
                line = samout.readline()
                continue

            if (flag & 0x40) > 0:
                if (flag & 0x10) == 0:
                    revkey = 'first-fwd'
                else:
                    revkey = 'first-rev'
            else:
                if (flag & 0x10) == 0:
                    revkey = 'second-fwd'
                else:
                    revkey = 'second-rev'

            parsedcigar = cigarparser.findall(cigar) # [(num, operator), (num, operator) ...], operator: [MDISX=]
            parsedcigar = [ (int(a), b) for a, b in parsedcigar ] # conver numbers in CIGAR to integers 
            
            #calculate rightmost genomic loc marked by cigar op
            #calculate rightmost read loc marked by cigar op
            #genomepos_rightmost = []
            #readpos_rightmost = []
            curgenomepos = pos - 1
            curreadpos = 0
            for num, op in parsedcigar:
                if op == 'N' or op == 'D': #not in read but reference
                    curgenomepos = curgenomepos + num
                    #genomepos_rightmost.append( curgenomepos )
                    #readpos_rightmost.append( curreadpos )
                elif op == 'I' or op == 'S': #not in reference but read
                    #genomepos_rightmost.append( curgenomepos+0.5 )
                    curreadpos = curreadpos + num
                    #readpos_rightmost.append( curreadpos )
                else: #M, X, =
                    curgenomepos = curgenomepos + num
                    #genomepos_rightmost.append( curgenomepos )
                    curreadpos = curreadpos + num
                    #readpos_rightmost.append( curreadpos )
                if curgenomepos >= varpos: #got it
                    if op == 'N' or op == 'D': #the variant is not in this read
                        numread = numread - 1
                        break
                    varreadpos = curreadpos - curgenomepos + varpos 
                    if varreadpos-1 < 0 or varreadpos > len(seq):
                        print varchrm, varpos, varreadpos, flag, pos, cigar, len(seq), curgenomepos, curreadpos
                        exit()
                    if seq[varreadpos - 1] == alt:
                        counter[revkey][varreadpos] += 1
                    break

            line = samout.readline()
        #done counting for this variant
        closefiles([samout, samerr])
        outputbuffer = outputbuffer + outputline( varchrm, varpos, numread, counter, readlength ) 
        if nvariant%500 == 0:
            out.write(outputbuffer)
            outputbuffer = ''
    if outputbuffer != '':
        out.write(outputbuffer)
    closefiles([vcf, out])
except:
    closefiles([vcf, samout, samerr, out])
    traceback.print_exc()

print time.time()-start
