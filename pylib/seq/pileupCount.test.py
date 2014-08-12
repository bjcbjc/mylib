#!/usr/bin/python

import argparse
import subprocess
import traceback
import re
from collections import Counter, OrderedDict, namedtuple
from randstr import randstr
from baseStrLib import BaseString
import time

def chrmReplaceFunc(matchobj):
    if matchobj.group('chrm') == 'MT': return 'chrM'
    elif matchobj.group('chrm') == 'chrM': return 'MT'
    elif 'chr' not in matchobj.group('chrm'): return 'chr'+matchobj.group('chrm')
    else: return matchobj.group('chrm')[3:]
    

    
def createTempLocFile(args, validChrms):
    #if vcf file, create a temporary file
    f = open(args.f)
    chrmformat = f.readline().strip('\n>').split()[0]
    f.close()

    if chrmformat not in validChrms:
        print 'inconsistent chrm format, reference: %s, bam: %s'%(chrmformat, validChrms[0])
        exit(1)
        
    f = open(args.loc)
    line = f.readline()
    f.close()
    removelater = []
    if 'fileformat=VCF' in line:
        locfile = args.tmpdir + randstr() + '.loc'
        #create temporary location file
        subprocess.call('grep -v ^# %s | cut -f1,2 > %s'%(args.loc, locfile), shell=True)
        removelater.append( locfile )
    elif len(line.split('\t')) > 2:
        locfile = args.tmpdir + randstr() + '.loc'
        #create temporary location file
        subprocess.call('cut -f1,2 %s > %s'%(args.loc, locfile), shell=True)
        removelater.append( locfile )
    else:
        locfile = args.loc

    #check location chrm format and filter unavailable contigs
    f = open(locfile)
    loc = f.readlines()
    f.close()
    locChrms = list(set(map(lambda(l): l.split()[0], loc)))    
    
    bamWithChr, locWithChr = False, False
    reprog = re.compile('^chr')
    if any(map(lambda(l): reprog.match(l), validChrms)):
        bamWithChr = True
        replaceProg = re.compile('(?P<chrm>^[0-9]{1,2}|^[XYMT]{1,2})')
    if any(map(lambda(l): reprog.match(l), locChrms)):
        locWithChr = True
        replaceProg = re.compile('(?P<chrm>^chr)')
    if bamWithChr != locWithChr:
        loc = [ re.sub(replaceProg, chrmReplaceFunc, line) for line in loc]

    loc = filter(lambda(line): line.split()[0] in validChrms, loc)
    #write another temporary file
    locfile = args.tmpdir + randstr() + '.loc'
    print 'tmpfile: %d lines, '%len(loc), locfile
    f = open(locfile, 'w')
    for line in loc: f.write(line)
    f.close()
    removelater.append( locfile )
    
    return locfile, removelater

def manualFilterBase(readdata, args, para):
    #one read
    #readdata: chr, pos, ref, #read, read_bases, read_qual, pos_in_read (if ignoreEnd is specified)
    #readdata[4] = para.baseCleanProg.sub('', readdata[4]) 
    readdata[4] = BaseString.splitBaseString(readdata[4]) #need this because of indels

    test = []
    if para.manualFilterBaseQ > 0:
        test.append( [ (ord(b) - para.baseQualOffset) >= para.manualFilterBaseQ for b in readdata[5] ] )
    if args.ignoreEnd > 0:
        test.append( [ min(args.readLength-p, p) > args.ignoreEnd for p in readdata[6] ] )

    passed = test[0]
    for i in range(1, len(test)):
        passed = [ a and b for a,b in zip(passed, test[i]) ]
    for i in range(4,len(readdata)):
        if i == 6:
            readdata[i] = [ b for b,p in zip(readdata[i],passed) if p ] 
        else:
            readdata[i] = ''.join( [ b for b,p in zip(readdata[i],passed) if p ] )
    readdata[3]  = '%d'%len(readdata[5]) #don't use readdata[4] beacause it may contain indels
    return readdata

def pileupToCount(samout, out, args, para, singlesite=False):
    #pileup special chars in read_bases: ^: start of read, followed by an additional char for mapping qual
    # $: end of read
    # > or < suggests non-coding region
    excludePattern = re.compile('(\^.)|(\$)')

    if args.tableformat and not singlesite:
        out.write('\t'.join(['chr', 'pos', 'ref', '#read']) + '\t' + '\t'.join(para.baseLabels) + '\n')

    line = samout.readline()
    count = 1
    outbuffer = ''
    while line:
        line = line.split() #chr, pos, ref, #read, read_bases, read_qual, pos_in_read (if ignoreEnd is specified)

        if int(line[3]) < args.r: 
            line = samout.readline()
            continue
        line[4] = excludePattern.sub('', line[4])

        if args.ignoreEnd > 0:
            line[6] = [ int(p) for p in line[6].split(',') ]
        if para.manualFilterBaseQ > 0 or args.ignoreEnd > 0:
            line = manualFilterBase(line, args, para)

        if int(line[3]) >= args.r: #need to check again because reads might be filtered out if manualFilterBaseQ is turned on
            if not args.stranded: 
                line[4] = line[4].upper().replace(',','.').replace('<','>')
            
            line[4] = line[4].replace('.', line[2].upper()).replace(',', line[2].lower())
            indels, line[4] = BaseString.matchIndels(line[4])

            #count
            #baseCount = Counter( re.sub(excludePattern, '', line[4] ))   #dict with counts and letters
            baseCount = Counter( line[4] )   #dict with counts and letters
            baseCount.update(indels)

            if args.tableformat:
                countstring = '\t'.join( [ '%s'%baseCount[b] for b in para.baseLabels if 'indel' not in b ] )
                if args.stranded:
                    indelKey = set( [ k.upper() for k in indels.keys() ] )
                    if args.indelMinReadBothStrand > 0:
                        indelKey = [ k for k in indelKey if min(indels[k], indels[k.lower()]) >= args.indelMinReadBothStrand ]
                    if len(indelKey) > 0:
                        indelCountString = '\t'.join( [ ','.join(indelKey), ','.join( ['%s'%indels[k] for k in indelKey] ), ','.join( ['%s'%indels[k.lower()] for k in indelKey] ) ] )
                    else:
                        indelCountString = '0\t0\t0'
                else:
                    indelKey = indels.keys()
                    if len(indelKey) > 0:
                        indelCountString = '\t'.join( [ ','.join(indelKey), ','.join( ['%s'%indels[k] for k in indelKey] ) ] )
                    else:
                        indelCountString = '0\t0'
                #out.write('\t'.join(line[:4]) + '\t' + countstring + '\n')
                outbuffer = outbuffer + '\t'.join(line[:4]) + '\t' + countstring + '\t' + indelCountString + '\n'
            else:
                baseCount = OrderedDict( sorted(baseCount.items(), key=lambda t:t[0]) )
                #decide what to write to the output: chr, pos, ref, total_read, base_count_string
                baseCountString = ''.join([ k+'%s'%v for k, v in baseCount.iteritems()])
                #out.write('\t'.join(line[:4]) + '\t' + baseCountString + '\n')            
                outbuffer = outbuffer + '\t'.join(line[:4]) + '\t' + baseCountString + '\n'

            count = count + 1

        if count%100 == 0:
            out.write(outbuffer)
            outbuffer = ''
        line = samout.readline()
    if outbuffer != '' and not singlesite:
        out.write(outbuffer)
        outbuffer = ''
    return outbuffer

def bamChrms(args):
    sampipe = subprocess.Popen(args.samtools + ' view -H ' + args.bam, shell=True, stdout=subprocess.PIPE)
    samheader = sampipe.stdout.read().split('\n')
    sampipe.stdout.close()
    samheader = filter(lambda(l): l[:3]=='@SQ', samheader)
    chrms = map(lambda(l): l.split()[1].replace('SN:',''), samheader)
    return chrms

def closefiles(files):
    for f in files:
        if f != None: f.close()
    return

argp = argparse.ArgumentParser(prog='pileupCount.py')
argp.add_argument('-d', default=10000, metavar='INT', type=int, help='INT for samtools mpileup, read maximally INT reads per bam [10,000]' )
argp.add_argument('-q', default=0, metavar='INT',  type=int, help='INT for samtools mpileup, min mapping quality for an alignment to be used [0]' )
argp.add_argument('-Q', default=13, metavar='INT', type=int, help='INT for samtools mpileup, mim base quality for a base to be considered [13]. NOTE: By default, samtools calculates BAQ (base alignment quality) and this threshold is applied to BAQ instead of raw base quality. To apply this threshold to raw base quality, specify -B to disable samtools from calculating BAQ.' )
argp.add_argument('-f', type=str, default='/nethome/bjchen/Projects/Simon/h37_hg19_chrOnly.fa', metavar='file', help='faidx-indexed reference fasta' )
argp.add_argument('-loc', type=str, metavar='file', default='', help='file, file can be VCF or bed, [chr, pos] as columns')
argp.add_argument('-reg', type=str, metavar='file', default='',help='file, file contains regions for allele counts; mutual exclusive with -loc; one output file will be generated by each region')
argp.add_argument('-bam', type=str, metavar='file', required=True, help='file, bam file for pileup')
argp.add_argument('-r',  metavar='INT', required=True, type=int, default=10, help='min reads required for each base')
argp.add_argument('-o', type=str, default='test.out', metavar='file', help='file name for the output, [<bam>.out]')
argp.add_argument('-tableformat', action='store_true', help='make table output instead; this option is turned off by default')
argp.add_argument('-stranded', action='store_true', help='a switch to specify if it is strand-specific data, [false]' )
argp.add_argument('-useSingleSiteCall', action='store_true', help='a switch to specify calling samtools with -r; this may be faster for sites with very high coverage, [false]' )
argp.add_argument('-ignoreEnd', default=0, metavar='INT', type=int, help='ignore INT bp from both end of each read [0]' )
argp.add_argument('-readLength', default=0, metavar='INT', type=int, help='used when -ignoreEnd is specified; if not -readLength is not specified, it will be inferred from a random read from the bam file [0]')
argp.add_argument('-Phred64', action='store_true', help='a switch to specify base quality encoding is Phred+64; default is Phred+33; [false]')
argp.add_argument('-indelMinReadBothStrand', default=0, metavar='INT', type=int, help='for table format output only; requires minimum number of reads from both forward and reverse strands in order for indels to be counted; this option is ignored if -stranded is not turned on [0]' )
#argp.add_argument('-ignore5', default=0, metavar='INT', type=int, help='ignore INT bp from 5'' end of each read')
#argp.add_argument('-ignore3', default=0, metavar='INT', type=int, help='ignore INT bp from 3'' end of each read')
argp.add_argument('-samtools', type=str, default='/data/NYGC/Software/samtools/samtools-0.1.19/samtools', help='samtools path, [/data/NYGC/Software/samtools/samtools-0.1.19]')
argp.add_argument('-tmpdir', type=str, default='./', help='path for temporary files, [./]')
argp.add_argument('-keepTmpLoc', action='store_true', help='keep temporary location file, for debug, [false]' )
#argp.add_argument('-refbasecheck', nargs='?', action='store_true', help='a switch to turn on double checking of reference bases in vcf and in genome; only valid if location file is in vcf format, [false]')
argp.add_argument('-samtools_mpileup_option', nargs='?', metavar='value', help='options that are passed to samtools mpileup')

args = argp.parse_known_args()
passon = ''
if isinstance(args, tuple):
    passon = ' '.join(args[1])
    args = args[0]
if args.tmpdir[-1] != '/': args.tmpdir = args.tmpdir + '/'

argtb = vars(args)
sampara = ['d', 'q', 'Q', 'f']

if args.loc != '' and args.reg != '':
    print 'Both -loc and -reg specified; ignore -loc.'
    args.loc = ''
elif args.loc == '' and args.reg == '':
    print 'Must specify -loc or -reg'
    exit()


PARA = namedtuple('PARA','baseLabels, manualFilterBaseQ, baseQualOffset')
para = PARA(baseLabels=['>', 'A', 'T', 'C', 'G', 'N', '*', 'indel', 'indelCount'], manualFilterBaseQ=0, baseQualOffset=33)
if args.stranded:
    para = para._replace(baseLabels = ['>', 'A', 'T', 'C', 'G', 'N', '<', 'a', 't', 'c', 'g', 'n', '*', 'indel', 'indelForwardCount','indelReverseCount'])
if args.ignoreEnd > 0:
    para = para._replace(manualFilterBaseQ = args.Q )
    args.Q = 0
    if '-O' not in passon: passon = passon + ' -O '
    #if '-B' not in passon: passon = passon + ' -B '
    if args.readLength == 0: 
        sampipe = subprocess.Popen(args.samtools + ' view ' + args.bam, shell=True, stdout=subprocess.PIPE)
        line = sampipe.stdout.readline().split('\t')
        sampipe.stdout.close()
        args.readLength = len(line[9])
        #print sampipe.pid
        sampipe.kill()
    args.readLength = args.readLength + 1
if args.Phred64:
    para = para._replace(baseQualOffset = 64)

validChrms = bamChrms(args)
if args.reg != '':
    regions = [l.strip() for l in open(args.reg).readlines()]
    if any([ l.split(':')[0] not in validChrms for l in regions ]):
        print 'Some chromosomes specified in regions are not found in bam'
        exit()
    cmd = args.samtools + ' mpileup -r %s '
    tmpfiles = []
else:
    locfile, tmpfiles = createTempLocFile(args, validChrms)
    if args.useSingleSiteCall:
        cmd = args.samtools + ' mpileup -r %s:%s-%s '
    else:
        cmd = args.samtools + ' mpileup -l %s '%(locfile)

for switch in sampara:
    cmd = cmd + '-%s %s '%(switch, argtb[switch])

cmd = cmd + ' ' + passon + ' ' + args.bam    
print 'Run command: ' + cmd

samout, samerr, out = None, None, None
try:
    #run samtools and pipe the output for further process (save space)
    if args.loc != '':
        if not args.useSingleSiteCall:
            sampipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            samout, samerr = sampipe.stdout, sampipe.stderr
            out = open(args.o, 'w')
            pileupToCount(samout, out, args, para, singlesite=False)
        else:
            siteloc = [l.strip().split() for l in open(locfile).readlines()]
            out = open(args.o, 'w')
            if args.tableformat:
                out.write('\t'.join(['chr', 'pos', 'ref', '#read']) + '\t' + '\t'.join(para.baseLabels) + '\n')
            count = 0
            outbuffer = ''
            for chrm, site in siteloc:
                #tic = time.time()
                sampipe = subprocess.Popen(cmd%(chrm,site,site), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                samout, samerr = sampipe.stdout, sampipe.stderr
            
                #out = open(args.o, 'w')
                outbuffer = outbuffer + pileupToCount(samout, out, args, para, singlesite=True)
                logmsg = samerr.read().strip()
                if logmsg != '[mpileup] 1 samples in 1 input files':
                    print logmsg
                count = count + 1
                #print chrm,site
                if count%100 == 0:
                    out.write(outbuffer)
                    outbuffer = ''
                closefiles([samout, samerr])
                #print '%d, %s'%(count, time.time()-tic)
            if outbuffer != '':
                out.write(outbuffer)
                outbuffer = ''
        if len(tmpfiles) > 0 and not args.keepTmpLoc:
            #remove temp file
            subprocess.call( 'rm -f ' + ' '.join(tmpfiles), shell=True )
        closefiles([samout, samerr, out])
    else: #regions; call mpileup with one region at a time; this is faster than -l
        for subr in regions:
            sampipe = subprocess.Popen(cmd%(subr), shell=True, stdout=subprocess.PIPE)#, stderr=subprocess.PIPE)
            samout = sampipe.stdout #, sampipe.stderr

            subr = subr.replace(',','')
            outfntag = '.'.join( re.split('[:-]', subr) )
            out = open(args.o + '.' + outfntag, 'w')
            pileupToCount(samout, out, args, para)
            closefiles([samout, samerr, out])
except:
    print 'Error somewhere, closing streaming pipes'
    if len(tmpfiles) > 0 and not args.keepTmpLoc:
        #remove temp file
        subprocess.call(  'rm -f ' + ' '.join(tmpfiles), shell=True )
    closefiles([samout, samerr, out])
    traceback.print_exc()

