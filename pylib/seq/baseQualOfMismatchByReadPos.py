

import argparse
import subprocess
import traceback
import re
import time
from collections import Counter

# Given a bam file and reference sequences, record (1) the distribution of reads containing
# mismatches, insertions, deletions, or soft clipping in the alignment and (2) number of
# mismatches, insertions, deletions, or soft clipping at each base position of reads.
#
# For now, alignments with any bit in 0xf04 are excluded from counting. This means we only count
# the primary alignments with properly-mapped pairs.
#
# The script calls samtools to calculate MD tag. Combining both CIGAR and MD, the statistics
# are collected. (Note: in CIGAR, M can be either match or mismatches; therefore MD is needed.)
# 
# Required: samtools, gawk
#
# Input: mismatchCountByRead.py -f reference.fa -bam alignment.bam -o output.filename [-samtools path/to/samtools] [-n] [-r read_length]
#
# Output: first line is the number of total counted reads. The following lines record the counts
# of reads for each category in (1) and (2) described above. Specifically, M, I, D and S represents
# mismatch, insertion, deletion, and soft-clipping, respectively. 'fwd' and 'rev' represents the counts are for 
# forward reads and reverse reads, respectively. Finally, rM/rI/rD/rS records the number of reads with
# xx number of M, I, D, or S (#1 above), whereas M/I/D/S records the positional counts (#2 above).
# Each column represents the number of M/I/D/S (1) or the position in reads (2). Indexes of columns start from 0.
#
# For example, assume the read length is 5.
#
#    fwd-rM  100  10  5  0  0  0
#    This means for forward reads, 100 reads have 0 mismatch, 10 have 1 mismatch, 5 have 2 mismatches, and
#    none have more than 2 mismatches.
#
#    rev-D  0  0  10  7  5  0
#    This means for reverse reads, no reads have deletion at position 1, 10 reads have deleletion at position 2,
#    7 have deletion at position 3, 5 have deletion at position 4, and 0 have deletion at position 5. The first
#    0 is a place holder (always zero, should be ignored).
#    
#

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
argp = argparse.ArgumentParser(prog='baseQualOfMismatchByReadPos.py')
#argp.add_argument('-f', type=str, default='/nethome/bjchen/Projects/Simon/h37_hg19_chrOnly.fa', metavar='reference_seq', help='faidx-indexed reference fasta' )
argp.add_argument('-input', type=str, metavar='bam_file', required=True, help='input file')
argp.add_argument('-o', type=str, default='test.mismatch.stat', metavar='output_file', help='file name for the output, [<bam>.mismatch.stat]')
argp.add_argument('-n', help='count N as mismatch if specified', action='store_true')

args = argp.parse_known_args()
passon = ''
if isinstance(args, tuple):
    passon = ' '.join(args[1])
    args = args[0]
print args

# Command string; "samtools view -F" is used to filter reads and then pipe into "samtools calmd" which is piped to gawk to retain only flag-information, CIGAR, and MD. The printed information (three columns per alignment) is then piped to this python script for parsing and counting (See below). 

# filter reads first to reduce the number of warnings in log, if MD is recalculated
#cmd = args.samtools + ' view -uh -F %d -q %d '%(0xf04, args.q) + args.bam + ' | ' + args.samtools + ' calmd  - ' + args.f  + ''' | gawk '{OFS="\t"; if ( and($2,0x2)>0 ) { for(i=1;i<=NF;i++) {if ($i ~ /MD:Z:/) {sub(/MD:Z:/,"",$i); print $2,$6,$i,$11}}} }' '''


# container for statistics to be collected
counter3 = {'first-fwd':Counter(), 'first-rev':Counter(), 'second-fwd':Counter(), 'second-rev':Counter()}
counter13 = {'first-fwd':Counter(), 'first-rev':Counter(), 'second-fwd':Counter(), 'second-rev':Counter()}
counter33 = {'first-fwd':Counter(), 'first-rev':Counter(), 'second-fwd':Counter(), 'second-rev':Counter()}
counter3_null = {'first-fwd':Counter(), 'first-rev':Counter(), 'second-fwd':Counter(), 'second-rev':Counter()}
counter13_null = {'first-fwd':Counter(), 'first-rev':Counter(), 'second-fwd':Counter(), 'second-rev':Counter()}
counter33_null = {'first-fwd':Counter(), 'first-rev':Counter(), 'second-fwd':Counter(), 'second-rev':Counter()}

readloc = [3, 13, 33]
start = time.time()
try:
    # regular expression parser for CIGAR and MD
    cigarparser = re.compile('(\d+)([MIDXS\=])')
    mdparser = re.compile('(\d+)[A-Z]')
    md_N_parser = re.compile('(?P<test>(?:\d+N)+\d+)')

    input = open(args.input)
    line = input.readline()
    while line:
        flag, cigar, md, baseQ3, baseQ13, baseQ33 = line.split()
        md = md.upper()
        flag = int(flag)

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

        if not args.n: # don't count N in MD as mismatach
            strN = md_N_parser.findall(md)  # (\d+)N(\d+)N...
            #count N as match and add it to the near-by matches (which are numbers in MD)
            replacement = [ sum(map(int,l.split('N')))+l.count('N') for l in strN ]
            replacement = map(str, replacement)

            #strN may have duplicates, but since len(strN) should be small, loop over all
            #the duplicated ones besides the first  will do nothing
            for i in range(len(strN)): 
                md = md.replace(strN[i], replacement[i])            
            cM = sum([ md.count(a) for a in ['A','T','C','G'] ]) #number of mismatches, don't count N
        else:
            cM = sum([ md.count(a) for a in ['A','T','C','G','N'] ]) #number of mismatches, count N

        parsedcigar = cigarparser.findall(cigar) # [(num, operator), (num, operator) ...], operator: [MDISX=]
        parsedmd = mdparser.findall(md) # [num, num, ...], num of matches
        parsedcigar = [ (int(a), b) for a, b in parsedcigar ] # conver numbers in CIGAR to integers 

        #exclude D for insertion positions
        parsedcigar_noD = [ l for l in parsedcigar if l[1]!='D' ]
        cigar_num, cigar_op = zip(*parsedcigar_noD) # CIGAR pos needs to exclude deletion (D) in reference, since these bases are not in the read
        insertion_num = [ num for num, op in zip(cigar_num, cigar_op) if op == 'I' ] # number of inserted bases
        insertion_pos = [ pos for pos, op in zip(cumsum(list(cigar_num)), cigar_op) if op == 'I' ] #pos of last inserted base


        #do positional mismatch counter, but md string doesn't account for insertion; so need cigar to offset the positions        
        hasM3, hasM13, hasM33 = False, False, False
        if len(parsedmd) > 0:
            parsedmd = map(int, parsedmd)
            mismatch_pos = [ a+b for a,b in zip(cumsum(parsedmd), range(1,len(parsedmd)+1)) ]
            for mispos in mismatch_pos:
                for insnum, inspos in zip(insertion_num, insertion_pos):
                    if mispos > inspos:
                        mispos = mispos + insnum
                if mispos == 3:
                    counter3[revkey][ord(baseQ3)-33] += 1
                    hasM3 = True
                elif mispos == 13:
                    counter13[revkey][ord(baseQ13)-33] += 1
                    hasM13 = True
                elif mispos == 33:
                    counter33[revkey][ord(baseQ33)-33] += 1
                    hasM33 = True

            if not hasM3:
                counter3_null[revkey][ord(baseQ3)-33] += 1
            if not hasM13:
                counter13_null[revkey][ord(baseQ13)-33] += 1
            if not hasM33:
                counter33_null[revkey][ord(baseQ33)-33] += 1
        else:
            counter3_null[revkey][ord(baseQ3)-33] += 1
            counter13_null[revkey][ord(baseQ13)-33] += 1
            counter33_null[revkey][ord(baseQ33)-33] += 1

        line = input.readline()

    
    # done counting; write output
    out = open(args.o, 'w')
    allcounter = [counter3,counter13,counter33,counter3_null,counter13_null,counter33_null]
    countername = ['counter3', 'counter13', 'counter33', 'counter3_null', 'counter13_null', 'counter33_null']
    for ctb, cname in zip(allcounter, countername):
        for direction in ctb:
            if len(ctb[direction]) > 0:
                out.write('%s\t%s'%(cname, direction))
                for ss in range(42):
                    out.write('\t%d'%ctb[direction][ss])
                out.write('\n')
    closefiles([ out])
except:
    closefiles([ out])
    traceback.print_exc()

print time.time()-start
