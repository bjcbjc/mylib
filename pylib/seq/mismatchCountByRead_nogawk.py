

import argparse
import subprocess
import traceback
import re
import time

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
# Required: samtools
#
# Input: mismatchCountByRead.py -f reference.fa -bam alignment.bam -o output.filename [-samtools path/to/samtools]
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
argp = argparse.ArgumentParser(prog='mismatchCountByRead.py')
argp.add_argument('-f', type=str, default='/nethome/bjchen/Projects/Simon/h37_hg19_chrOnly.fa', metavar='reference_seq', help='faidx-indexed reference fasta' )
argp.add_argument('-bam', type=str, metavar='bam_file', required=True, help='bam file')
argp.add_argument('-o', type=str, default='test.mismatch.stat', metavar='output_file', help='file name for the output, [<bam>.mismatch.stat]')
argp.add_argument('-n', help='count N as mismatch if specified', action='store_true')
argp.add_argument('-r', type=int, metavar='read_length', default=0, help='max read length; if not provided, obtained from the first read')
argp.add_argument('-samtools', type=str, metavar='path_to_samtools',default='/data/NYGC/Software/samtools/samtools-0.1.19/samtools', help='samtools path, [/data/NYGC/Software/samtools/samtools-0.1.19]')

args = argp.parse_known_args()
passon = ''
if isinstance(args, tuple):
    passon = ' '.join(args[1])
    args = args[0]
print args
# Command string; "samtools calmd" is called, then piped to "samtools view" to print out alignment, which is piped to gawk to retain only flag-information, CIGAR, and MD. The printed information (three columns per alignment) is then piped to this python script for parsing and counting (See below).
#cmd = args.samtools + ' calmd ' + args.bam + ' ' + args.f + ''' | gawk '{OFS="\t"; if ( and($2,0xf04)==0 ) { sub(/MD:Z:/,"",$NF); print and($2,0x10),$6,$NF} }' '''

# filter reads first to reduce the number of warnings in log, if MD is recalculated
cmd = args.samtools + ' view -uh -F %d '%(0xf04) + args.bam + ' | ' + args.samtools + ' calmd  - ' + args.f 

print 'Run command: ' + cmd

# set up pipes
samout, samerr, out = None, None, None

# container for statistics to be collected
counter = {'fwd':{'D':[], 'I':[], 'M':[], 'S':[], 'rD':[], 'rI':[], 'rM':[], 'rS':[]}, 'rev':{'D':[], 'I':[], 'M':[], 'S':[], 'rD':[], 'rI':[], 'rM':[], 'rS':[] }}

start = time.time()
try:
    # decide read length first
    sampipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    samout = sampipe.stdout #, sampipe.stderr

    line = samout.readline()
    while line[0] == '@':
        line = samout.readline()

    # now we have the first read
    if args.r == 0:
        readlength = len(line.split()[9])
    else:
        readlength = args.r

    # initialized counters
    for k in counter.keys():
        for token in counter[k].keys():
            counter[k][token] = [0.0]*(readlength+1)

    # regular expression parser for CIGAR and MD
    cigarparser = re.compile('(\d+)([MIDXS\=])')
    mdparser = re.compile('(\d+)[A-Z]')
    md_N_parser = re.compile('(?P<test>(?:\d+N)+\d+)')

    numread = 1
    while line:
        md = line.split('MD:Z:')[1].split()[0].upper() # faster than re
        line = line.split()
        cigar = line[5]
        revCmpFlag = (int(line[1]) & 0x10) > 0

        if cigar == '*': #not available
            line = samout.readline()
            continue

        revkey = 'fwd'
        if revCmpFlag:
            revkey = 'rev'

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

        cigar_num, cigar_op = zip(*parsedcigar)        
        deletion_num = [ num for num, op in parsedcigar if op == 'D'] # number of deleted bases (as in reference)
        deletion_pos = [ pos for pos, op in zip(cumsum(list(cigar_num)),cigar_op) if op == 'D' ]
        deletion_pos = [ pos-num for pos, num in zip(deletion_pos, cumsum(deletion_num)) ] # position of deletion

        #read counters
        for token in ['D', 'I', 'S']: #count deletion, insertion, soft-clipping
            c = sum(  [ num for num,op in parsedcigar if op==token ]  )
            counter[revkey]['r%s'%token][c] = counter[revkey]['r%s'%token][c] + 1
            if token == 'D':
                cM = cM - c #minus deletion
        counter[revkey]['rM'][cM] = counter[revkey]['rM'][cM] + 1

        #positional insertion counter
        for insnum, inspos in zip(insertion_num, insertion_pos):
            counter[revkey]['I'][inspos-insnum+1:inspos+1] = [ cur+1 for cur in counter[revkey]['I'][inspos-insnum+1:inspos+1] ]

        #positional deletion counter
        for delnum, delpos in zip(deletion_num, deletion_pos):
            w = delnum/2.0
            counter[revkey]['D'][delpos:delpos+2 ] = [ cur+w for cur in counter[revkey]['D'][delpos:delpos+2] ]

        #positional clipping counter
        if parsedcigar[0][1] == 'S':
            counter[revkey]['S'][1:parsedcigar[0][0]+1] = [ cur+1 for cur in counter[revkey]['S'][1:parsedcigar[0][0]+1] ]
        if parsedcigar[-1][1] == 'S':
            counter[revkey]['S'][readlength+1-parsedcigar[-1][0]:readlength+1] = [ cur+1 for cur in counter[revkey]['S'][readlength+1-parsedcigar[-1][0]:readlength+1] ]

        #do positional mismatch counter, but md string doesn't account for insertion; so need cigar to offset the positions        
        if len(parsedmd) > 0:
            parsedmd = map(int, parsedmd)
            mismatch_pos = [ a+b for a,b in zip(cumsum(parsedmd), range(1,len(parsedmd)+1)) ]
            for mispos in mismatch_pos:
                for insnum, inspos in zip(insertion_num, insertion_pos):
                    if mispos > inspos:
                        mispos = mispos + insnum
                counter[revkey]['M'][mispos] = counter[revkey]['M'][mispos] + 1

        line = samout.readline()
        numread = numread + 1
    
    # done counting; write output
    out = open(args.o, 'w')
    out.write('%d\n'%numread)
    for direction in counter:
        for optype in counter[direction]:
            out.write('%s-%s\t'%(direction, optype))
            out.write('\t'.join(map(str,counter[direction][optype])) + '\n')
    closefiles([samout, samerr, out])
except:
    closefiles([samout, samerr, out])
    traceback.print_exc()

print time.time()-start
