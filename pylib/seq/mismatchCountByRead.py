

import argparse
import subprocess
import traceback
import re

def cumsum(l):
    s = l[:] #make a copy
    for i in xrange(1, len(s)):
        s[i] = s[i]+s[i-1]
    return s

def closefiles(files):
    for f in files:
        if f != None: f.close()
    return

argp = argparse.ArgumentParser(prog='mismatchCountByRead.py')
argp.add_argument('-f', nargs='?', default='/nethome/bjchen/Projects/Simon/h37_hg19_chrOnly.fa', metavar='file', help='faidx-indexed reference fasta' )
argp.add_argument('-bam', nargs='?', metavar='file', required=True, help='bam file')
argp.add_argument('-o', nargs='?', default='test.mismatch.stat', metavar='file', help='file name for the output, [<bam>.mismatch.stat]')
argp.add_argument('-samtools', nargs='?', default='/data/NYGC/Software/samtools/samtools-0.1.19/samtools', help='samtools path, [/data/NYGC/Software/samtools/samtools-0.1.19]')

args = argp.parse_known_args()
passon = ''
if isinstance(args, tuple):
    passon = ' '.join(args[1])
    args = args[0]


cmd = args.samtools + ' calmd -u ' + args.bam + ' ' + args.f + ' | ' + args.samtools + ''' view - | gawk '{OFS="\t"; if ( and($2,0xf04)==0 ) { sub(/MD:Z:/,"",$NF); print and($2,0x10),$6,$NF} }' '''

print 'Run command: ' + cmd

samout, samerr, out = None, None, None
counter = {'fwd':{'D':[], 'I':[], 'M':[], 'S':[], 'rD':[], 'rI':[], 'rM':[], 'rS':[]}, 'rev':{'D':[], 'I':[], 'M':[], 'S':[], 'rD':[], 'rI':[], 'rM':[], 'rS':[] }}
try:
    #decide read length first
    sampipe = subprocess.Popen(args.samtools + ' view ' + args.bam, shell=True, stdout=subprocess.PIPE)
    seq = sampipe.stdout.readline()
    sampipe.stdout.close()
    readlength = len(seq.split('\t')[9])

    for k in counter.keys():
        for token in counter[k].keys():
            counter[k][token] = [0.0]*(readlength+1)

    #run samtools and pipe the output for further process (save space)
    sampipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)#, stderr=subprocess.PIPE)
    samout = sampipe.stdout #, sampipe.stderr

    cigarparser = re.compile('(\d+)([MIDXS\=])')
    mdparser = re.compile('(\d+)[A-Z]')
    line = samout.readline()
    numread = 1
    while line:
        revCmpFlag, cigar, md = line.split() #and(rev-cmp-flag), CIGAR, MD
        md = md.upper()
        revkey = 'fwd'
        if revCmpFlag != '0':
            revkey = 'rev'

        cM = sum([ md.count(a) for a in ['A','T','C','G','N'] ]) #number of mismatches
        parsedcigar = cigarparser.findall(cigar) # [(num, op), (num, op) ...], op:MDISX=
        parsedmd = mdparser.findall(md) # [num, num, ...], num of matches
        parsedcigar = [ (int(a), b) for a, b in parsedcigar ]

        #exclude D for insertion positions
        parsedcigar_noD = [ l for l in parsedcigar if l[1]!='D' ]
        cigar_num, cigar_op = zip(*parsedcigar_noD) #cigar pos needs to exclude deletion (D) in reference
        insertion_num = [ num for num, op in zip(cigar_num, cigar_op) if op == 'I' ]
        insertion_pos = [ pos for pos, op in zip(cumsum(list(cigar_num)), cigar_op) if op == 'I' ] #pos of last inserted base

        cigar_num, cigar_op = zip(*parsedcigar)        
        deletion_num = [ num for num, op in parsedcigar if op == 'D']
        deletion_pos = [ pos for pos, op in zip(cumsum(list(cigar_num)),cigar_op) if op == 'D' ]
        deletion_pos = [ pos-num for pos, num in zip(deletion_pos, cumsum(deletion_num)) ]

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

