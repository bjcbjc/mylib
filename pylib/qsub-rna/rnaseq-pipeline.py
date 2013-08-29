
from sys import argv
from os import system
import os.path
import copy
import cmdGenerator
import jobFactory
import configRobot




#read all available parameters for programs
def getAvailableParas():
    t = open('AllAvailableParas.txt').readlines()
    t = map(lambda(l):l.strip(), t)
    t = map(lambda(l):l.split('\t'), t)
    tb = {}
    for l in t:
        tb[ l[0] ] = l[1:]
    return tb

def DESeqPair(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmdset = configRobot.makeParasList(cmdset, ['meta', 'group1', 'group2'])
    cmd, mem, time, prefix = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'prefix'])
    group1, group2 = configRobot.popParas(cmdset, ['group1', 'group2'])
    template = open(cmdset.pop('template')).read()
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    cmdset['inputpath'] = inputpath
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    cmdset['outputpath'] = outputpath
    cmdset['prefix'] = prefix

    meta = configRobot.popParas(cmdset, ['meta'])
    if meta[0] == "''" and len(meta) == 1:
        cmdset['meta'] = 'c()'
    else:
        cmdset['meta'] = 'c(\'' + '\', \''.join(meta) + '\')'

    if cmdset['countfnprefix'] == "''": cmdset['countfnprefix'] = ''
    if cmdset['countfnsuffix'] == "''": cmdset['countfnsuffix'] = ''

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    setuppathcmd = cmdGenerator.formatCmd('source ~/libraries/setup_seqtools')
    for i in range(len(group1)):
        paraset = copy.deepcopy(cmdset)        
        paraset['gr1'] = group1[i]
        paraset['gr2'] = group2[i]
        jobfnprefix = prefix + '_' + group1[i] + '_' + group2[i]

        f = open('./%s.R'%(jobfnprefix), 'w')
        f.write(template%paraset)
        f.close()

        deseqcmd = cmdGenerator.formatCmd('Rscript', './%s.R'%(jobfnprefix))
        mvRscriptcmd = cmdGenerator.formatCmd('mv ./%s.R %s'%(jobfnprefix, outputpath))
        mvscriptcmd = cmdGenerator.formatCmd('mv ./%s%s %s'%(jobfnprefix, jobmanager.ext, outputpath))

        jobmanager.createJob(jobfnprefix, [setuppathcmd, deseqcmd, mvRscriptcmd, mvscriptcmd], outpath = outputpath, outfn = jobfnprefix)
    return jobmanager
    
    

def HTSeqCount(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, bam = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'bam'])
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath) 
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    gtf = configRobot.popParas(cmdset, ['GTF'])

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    setuppathcmd = cmdGenerator.formatCmd('source ~/libraries/setup_seqtools')

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]

    for sample in samples:
        paraset = copy.deepcopy(cmdset)        
        jobfnprefix = prefix + '_' + sample
        if bam == '=sample':
            inputfile = sample
            outputfile = sample + '.count'
        else: 
            inputfile = sample + '/' + bam
            outputfile = sample + '/' + bam + '.count'
            cmdGenerator.checkPath(outputpath + sample, create=createpath)

        samcmd = 'samtools view -h %s | '%(inputpath+inputfile)
        htseq = 'python -m HTSeq.scripts.count -q '
        countcmd = cmdGenerator.formatCmd(samcmd, htseq, paraset, '-', gtf, ' > %s'%(outputpath + outputfile))
        mvscriptcmd = cmdGenerator.formatCmd('mv ./%s%s %s'%(jobfnprefix, jobmanager.ext, outputpath))

        jobmanager.createJob(jobfnprefix, [setuppathcmd, countcmd, mvscriptcmd], outpath = outputpath, outfn = jobfnprefix)
    return jobmanager
        

#mapping
#pipeline component: run tophat to map data
def tophat(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix'])
    paired = configRobot.popParas(cmdset, ['paired'])
    readext, outpath, genome = configRobot.popParas(cmdset, ['readext', 'outputpath', 'genome'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    setuppathcmd = cmdGenerator.formatCmd('source ~/libraries/setup_seqtools\necho $BOWTIE2_INDEXES')

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]

    cmdGenerator.checkPath(outpath + '%s/'%prefix, create=createpath)

    for sample in samples:
        if paired == 'paired' or paired == 'yes':
            reads = map(lambda(i): inputpath + sample + '_%d'%i + readext, [1,2])
        elif paired == 'single' or paired == 'no':
            reads = inputpath + sample + readext
        paraset = copy.deepcopy(cmdset)
        paraset['-o'] = outpath + '%s/'%prefix + sample
        jobfnprefix = prefix + '_' + sample
        
        tophatcmd = cmdGenerator.formatCmd(cmd, paraset, genome, reads)
        mvscriptcmd = cmdGenerator.formatCmd('mv ./%s%s %s'%(jobfnprefix, jobmanager.ext, paraset['-o']))

        if int(paraset['-p']) > 1: #multiple threads per job
            sgeopt = ['-pe smp ' + paraset['-p']]
        else:
            sgeopt = []

        #need to create the output directory first, otherwise SGE complains cannot put the stdout in its path
        cmdGenerator.checkPath(paraset['-o'], create=createpath)
        jobmanager.createJob(jobfnprefix, [setuppathcmd, tophatcmd, mvscriptcmd], outpath = paraset['-o'], outfn = jobfnprefix, sgeopt=sgeopt)
    return jobmanager


def cufflinks(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, bam = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'bam'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    outputpath = cmdGenerator.checkPath(outputpath + '%s/'%prefix, create=createpath)

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]
    
    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    for sample in samples:
        jobname = prefix + '_' + sample
        CMD = []
        CMD.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )

        paraset = copy.deepcopy(cmdset)
        paraset['-o'] = outputpath + sample
        paraset = configRobot.validParas(paraset, availParas['cufflinks'])
        cmdGenerator.checkPath(paraset['-o'], create=createpath)
        CMD.append( cmdGenerator.formatCmd(cmd, paraset, inputpath+'%s/'%sample+bam) )
        
        CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, paraset['-o'])) )
        sgeopt = []
        if '-p' in paraset.keys():
            if int(paraset['-p']) > 1: #multi threads
                sgeopt = ['-pe smp ' + paraset['-p']]
        elif '--num-threads' in paraset.keys():
            if int(paraset['--num-threads']) > 1:
                sgeopt = ['-pe smp ' + paraset['-p']]
        jobmanager.createJob(jobname, CMD, outpath=paraset['-o'], outfn=jobname, sgeopt=sgeopt)
    return jobmanager


def cuffcompare(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, gtf = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'gtf'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]
    
    sampletext = ''
    for sample in samples:
        sampletext = sampletext + '%s%s/%s '%(inputpath, sample, gtf)

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    jobname = prefix
    CMD = []
    CMD.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )

    paraset = copy.deepcopy(cmdset)
    paraset['-o'] = outputpath + paraset['-o']
    paraset = configRobot.validParas(paraset, availParas['cuffcompare'])
    CMD.append( cmdGenerator.formatCmd(cmd, paraset, sampletext) )
    CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )

    jobmanager.createJob(jobname, CMD, outpath=outputpath, outfn=jobname, trackcmd=False)
    return jobmanager



def cuffmerge(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, gtf = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'gtf'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]
    
    sampletext = '"'
    for sample in samples:
        sampletext = sampletext + '%s%s/%s\\n'%(inputpath, sample, gtf)
    sampletext = sampletext + '"'

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    jobname = prefix
    CMD = []
    CMD.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )
    paraset = copy.deepcopy(cmdset)
    paraset['-o'] = outputpath + paraset['-o']

    CMD.append( cmdGenerator.formatCmd('echo', sampletext, '>', paraset['-o'] + '.samples') )
    
    paraset = configRobot.validParas(paraset, availParas['cuffmerge'])
    cmdGenerator.checkPath(paraset['-o'], create=createpath)
    CMD.append( cmdGenerator.formatCmd(cmd, paraset, paraset['-o'] + '.samples') )
    CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, paraset['-o'])) )
    CMD.append( cmdGenerator.formatCmd('rm -f', paraset['-o'] + '.samples') )
    sgeopt = []
    if '-p' in paraset.keys():
        if int(paraset['-p']) > 1: #multi threads
            sgeopt = ['-pe smp ' + paraset['-p']]
    elif '--num-threads' in paraset.keys():
        if int(paraset['--num-threads']) > 1:
            sgeopt = ['-pe smp ' + paraset['-p']]
    jobmanager.createJob(jobname, CMD, outpath=paraset['-o'], outfn=jobname, sgeopt=sgeopt, trackcmd=False)
    return jobmanager

def cuffdiff(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, gtf, bam = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'gtf', 'bam'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]
    
    sampletext = ''
    for sample in samples:
        sampletext = sampletext + '%s%s/%s '%(inputpath, sample, bam)

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    jobname = prefix
    CMD = []
    CMD.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )

    paraset = copy.deepcopy(cmdset)
    paraset = configRobot.validParas(paraset, availParas['cuffdiff'])
    cmdGenerator.checkPath(paraset['--output-dir'], create=createpath)
    CMD.append( cmdGenerator.formatCmd(cmd, paraset, gtf, sampletext) )
    CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, paraset['--output-dir'])) )
    sgeopt = []
    if '-p' in paraset.keys():
        if int(paraset['-p']) > 1: #multi threads
            sgeopt = ['-pe smp ' + paraset['-p']]
    elif '--num-threads' in paraset.keys():
        if int(paraset['--num-threads']) > 1:
            sgeopt = ['-pe smp ' + paraset['-p']]
    jobmanager.createJob(jobname, CMD, outpath=paraset['--output-dir'], outfn=jobname, sgeopt=sgeopt, trackcmd=False)
    return jobmanager


def cuffdiff_v1(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, gtf, bam = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'gtf', 'bam'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]
    
    sampletext = ''
    for sample in samples:
        sampletext = sampletext + '%s%s/%s '%(inputpath, sample, bam)

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    jobname = prefix
    CMD = []
    CMD.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )

    paraset = copy.deepcopy(cmdset)
    paraset = configRobot.validParas(paraset, availParas['cuffdiff'])
    cmdGenerator.checkPath(paraset['--output-dir'], create=createpath)
    CMD.append( cmdGenerator.formatCmd('/ifs/home/c2b2/dp_lab/bc2252/SeqTool/cufflinks_1/cuffdiff', paraset, gtf, sampletext) )
    CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, paraset['--output-dir'])) )
    sgeopt = []
    if '-p' in paraset.keys():
        if int(paraset['-p']) > 1: #multi threads
            sgeopt = ['-pe smp ' + paraset['-p']]
    elif '--num-threads' in paraset.keys():
        if int(paraset['--num-threads']) > 1:
            sgeopt = ['-pe smp ' + paraset['-p']]
    jobmanager.createJob(jobname, CMD, outpath=paraset['--output-dir'], outfn=jobname, sgeopt=sgeopt, trackcmd=False)
    return jobmanager



def countmismatches(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True
    
    cmd, mem, time, samples, prefix, bam = configRobot.popParas(cmdset,['cmd', 'mem', 'time', 'sample', 'prefix', 'bam'])
    
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    direct = ['forward', 'reverse']
    for sample in samples:
        for d in direct:
            jobname = prefix+'_'+sample+'_'+d
            CMD = []
            CMD.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )
            CMD.append( cmdGenerator.formatCmd('python','countmismatches.py',inputpath+sample+'/'+bam, outputpath+sample+'.%s.countmis'%d, d) )
            CMD.append( cmdGenerator.formatCmd('mv', '%s%s'%(jobname, jobmanager.ext), outputpath) )
            jobmanager.createJob(jobname, CMD, outpath=outputpath, outfn=jobname, trackcmd=False)
    return jobmanager
    

def filetersingleton(cmdset, runmode='test'): 
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, bam = configRobot.popParas(cmdset,['cmd', 'mem', 'time', 'sample', 'prefix', 'bam'])
    if 'TMP_DIR' in cmdset.keys():
        TMP_DIR = cmdset.pop('TMP_DIR')
    else:
        TMP_DIR = ''

    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    setuppathcmd = cmdGenerator.formatCmd('source ~/libraries/setup_seqtools')


    javacmd = 'java -Xmx%dg -jar'%(int(mem.replace('G',''))-1)
    samview = 'samtools view -b -h -F 8'
    reorder = 'ReorderSam.jar VALIDATION_STRINGENCY=LENIENT'
    RG = 'AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=LENIENT RGLB=dUTP RGPL=illumina RGPU=1'
    mdup = 'MarkDuplicates.jar'
    

    for sample in samples:
        jobfnprefix = prefix + '_' + sample
        paraset = copy.deepcopy(cmdset)
        if TMP_DIR != '': paraset['TMP_DIR'] = '%s'%TMP_DIR
        tmpbam = []

        CMDs = []
        CMDs.append(setuppathcmd)

        paraset['INPUT'] = '=%s/%s'%(inputpath+sample, bam)
        paraset['OUTPUT'] = '=%s/%s.reorder.bam'%(inputpath+sample, bam.replace('.bam',''))
        tmpbam.append(paraset['OUTPUT'].strip('='))
        CMDs.append( cmdGenerator.formatCmd(javacmd, programpath+reorder, paraset) )
        
        paraset['INPUT'] = '=%s/%s.reorder.bam'%(inputpath+sample, bam.replace('.bam',''))
        paraset['OUTPUT'] = '=%s/%s.reorder.addRG.bam'%(inputpath+sample, bam.replace('.bam',''))
        paraset['RGSM'] = '=%s'%sample        
        CMDs.append( cmdGenerator.formatCmd(javacmd, programpath+RG, paraset) )

        paraset = copy.deepcopy(cmdset)
        paraset['-o'] = '%s/%s.filter.bam'%(inputpath+sample, bam.replace('.bam','.addRG'))
        CMDs.append( cmdGenerator.formatCmd(samview, paraset, inputpath+sample+'/'+bam.replace('.bam','.addRG.bam')) )

        
        paraset = copy.deepcopy(cmdset)
        paraset['INPUT'] = '=%s/%s.filter.bam'%(inputpath+sample, bam.replace('.bam','.addRG'))
        paraset['OUTPUT'] = '%s/%s.mdup.bam'%(paraset['INPUT'].replace('.bam', ''))
        CMDs.append( cmdGenerator.formatCmd(javacmd, programpath+mdup, paraset) )

        paraset = copy.deepcopy(cmdset)
        CMDs.append( cmdGenerator.formatCmd('samtools index', bam.replace('.bam', '.reorder.addRG.filter.mdup.bam')) )


        CMDs.append( cmdGenerator.formatCmd('rm -f', tmpbam) )
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobfnprefix, jobmanager.ext, inputpath+sample)) )

        jobmanager.createJob(jobfnprefix, CMDs, outpath = inputpath+sample, outfn = jobfnprefix)
    return jobmanager


def markDup(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix = configRobot.popParas(cmdset,['cmd', 'mem', 'time', 'sample', 'prefix'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))
    bam = configRobot.popParas(cmdset, 'bam')
    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    javacmd = 'java -Xmx%dg -jar'%(int(mem.replace('G',''))-1)
    mdupjar = 'MarkDuplicates.jar'
    idxcmd = 'samtools index'

    for sample in samples:
        jobname = prefix + '_' + sample
        paraset = copy.deepcopy(cmdset)
        paraset['INPUT'] = '=%s/%s'%(inputpath + sample, bam)
        paraset['OUTPUT'] = paraset['INPUT'].replace('.bam', '.mdup.bam')
        paraset['METRICS_FILE'] = '=%s/%s'%(inputpath + sample, prefix + '_mdupmetrics.txt')
        paraset = configRobot.validParas(paraset, availParas[mdupjar])
        CMDs = []
        CMDs.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )
        CMDs.append( cmdGenerator.formatCmd(javacmd, programpath+mdupjar, paraset) )
        CMDs.append( cmdGenerator.formatCmd(idxcmd, paraset['OUTPUT'].strip('=')) )
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, inputpath+sample)) )
        jobmanager.createJob(jobname, CMDs, outpath = inputpath+sample, outfn = jobname)
    return jobmanager


def preGATK(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True
    
    cmd, mem, time, samples, prefix = configRobot.popParas(cmdset,['cmd', 'mem', 'time', 'sample', 'prefix'])
    
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    picardpath = cmdGenerator.checkPath(cmdset.pop('picardpath'))
    gatkpath = cmdGenerator.checkPath(cmdset.pop('gatkpath'))

    bam = configRobot.popParas(cmdset, 'bam')
    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    if '-Djava.io.tmpdir' in cmdset.keys():
        javacmd = 'java ' + '-Djava.io.tmpdir' + cmdGenerator.checkPath(cmdset.pop('-Djava.io.tmpdir'))
    else:
        javacmd = 'java'
    javacmd = javacmd + ' -Xmx%dg -jar'%(int(mem.replace('G',''))-2)
        

    samview = 'samtools view -b -h -F 264'
    reorder = picardpath + 'ReorderSam.jar VALIDATION_STRINGENCY=LENIENT'
    RG = picardpath + 'AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=LENIENT RGLB=dUTP RGPL=illumina RGPU=1'
    mdupjar = picardpath + 'MarkDuplicates.jar'
    GATK = gatkpath + 'GenomeAnalysisTK.jar '
    createTg = '-T RealignerTargetCreator '
    realign = '-T IndelRealigner '

    idxcmd = 'samtools index'
    clearup = 'rm -f '


    for sample in samples:
        CMDs = []
        CMDs.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )

        jobname = prefix + '_' + sample

        #filter
        paraset = copy.deepcopy(cmdset)
        paraset['-o'] = '%s/%s.filter.bam'%(inputpath+sample, bam.replace('.bam',''))
        lastoutput = paraset['-o']
        del paraset['-R']
        del paraset['-filterMBQ']
        #paraset = configRobot.validParas(paraset, availParas['samtools'])
        CMDs.append( cmdGenerator.formatCmd(samview, paraset, inputpath+sample+'/'+bam ) )

        #reorder by chrm
        paraset = copy.deepcopy(cmdset)
        paraset['INPUT'] = '=%s'%lastoutput
        paraset['OUTPUT'] = '=%s.reorder.bam'%(lastoutput.replace('.bam',''))
        paraset['REFERENCE'] = '=%s'%paraset['-R']
        paraset = configRobot.validParas(paraset, availParas['ReorderSam.jar'])
        CMDs.append( cmdGenerator.formatCmd(javacmd, reorder, paraset) )
        CMDs.append( cmdGenerator.formatCmd(clearup, lastoutput) )
        lastoutput = paraset['OUTPUT'].strip('=')
        
        #add RG
        paraset = copy.deepcopy(cmdset)
        paraset['INPUT'] = '=%s'%lastoutput
        paraset['OUTPUT'] = '=%s.addRG.bam'%(lastoutput.replace('.bam',''))
        paraset['RGSM'] = '=%s'%sample        
        paraset = configRobot.validParas(paraset, availParas['AddOrReplaceReadGroups.jar'])
        CMDs.append( cmdGenerator.formatCmd(javacmd, RG, paraset) )
        CMDs.append( cmdGenerator.formatCmd(clearup, lastoutput) )
        lastoutput = paraset['OUTPUT'].strip('=')

        #mark duplicates
        paraset = copy.deepcopy(cmdset)
        paraset['INPUT'] = '=%s'%lastoutput
        paraset['OUTPUT'] = '=%s.mdup.bam'%(lastoutput.replace('.bam', ''))
        paraset['METRICS_FILE'] = '=%s/%s'%(inputpath + sample, prefix + '_mdupmetrics.txt')
        paraset = configRobot.validParas(paraset, availParas['MarkDuplicates.jar'])
        CMDs.append( cmdGenerator.formatCmd(javacmd, mdupjar, paraset) )
        CMDs.append( cmdGenerator.formatCmd(idxcmd, paraset['OUTPUT'].strip('=')) )
        lastoutput = paraset['OUTPUT'].strip('=')

        #create intervals
        paraset = copy.deepcopy(cmdset)
        paraset['-I'] = lastoutput
        paraset['-o'] = lastoutput.replace('.bam', '.intervals')
        CMDs.append( cmdGenerator.formatCmd(javacmd, GATK+createTg, paraset) )

        #realign
        paraset['-targetIntervals'] = paraset['-o']
        paraset['-o'] = lastoutput.replace('.bam', '.realign.bam')
        CMDs.append( cmdGenerator.formatCmd(javacmd, GATK+realign, paraset) )

        #clear up
        CMDs.append( cmdGenerator.formatCmd(clearup, lastoutput) )
        CMDs.append( cmdGenerator.formatCmd(clearup, lastoutput.replace('.bam', '.intervals')) )
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, inputpath+sample)) )
        jobmanager.createJob(jobname, CMDs, outpath = inputpath+sample, outfn = jobname)
    
    return jobmanager


def picardReorderSam(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True
    
    cmd, mem, time, samples, prefix = configRobot.popParas(cmdset,['cmd', 'mem', 'time', 'sample', 'prefix'])
    
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'))
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))

    bam = configRobot.popParas(cmdset, 'bam')
    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    if '-Djava.io.tmpdir' in cmdset.keys():
        javacmd = 'java ' + '-Djava.io.tmpdir' + cmdGenerator.checkPath(cmdset.pop('-Djava.io.tmpdir'))
    else:
        javacmd = 'java'
    javacmd = javacmd + ' -Xmx%dg -jar'%(int(mem.replace('G',''))-2)        
    reorder = programpath + 'ReorderSam.jar VALIDATION_STRINGENCY=LENIENT'


    for sample in samples:
        CMDs = []
        CMDs.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )

        jobname = prefix + '_' + sample

        #reorder by chrm
        paraset = copy.deepcopy(cmdset)
        if bam == '=sample':
            inputfile = sample
        else:
            inputfile = sample + '/' + bam
        paraset['INPUT'] = '=%s'%(inputpath + inputfile)
        paraset['OUTPUT'] = '=%s.reorder.bam'%(outputpath + inputfile.replace('.bam',''))
        paraset = configRobot.validParas(paraset, availParas['ReorderSam.jar'])
        CMDs.append( cmdGenerator.formatCmd(javacmd, reorder, paraset) )
        
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )
        jobmanager.createJob(jobname, CMDs, outpath = outputpath, outfn = jobname)
    
    return jobmanager


def GATK_genotyper(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True
    
    cmd, mem, time, samples, prefix = configRobot.popParas(cmdset,['cmd', 'mem', 'time', 'sample', 'prefix'])
    
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    gatkpath = cmdGenerator.checkPath(cmdset.pop('gatkpath'))

    bam = configRobot.popParas(cmdset, 'bam')
    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    if '-Djava.io.tmpdir' in cmdset.keys():
        javacmd = 'java ' + '-Djava.io.tmpdir' + cmdGenerator.checkPath(cmdset.pop('-Djava.io.tmpdir'))
    else:
        javacmd = 'java'
    javacmd = javacmd + ' -Xmx%dg -jar'%(int(mem.replace('G',''))-2)
        
    GATK = gatkpath + 'GenomeAnalysisTK.jar '
    genotyper = '-T UnifiedGenotyper '


    CMDs = []
    CMDs.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )
    
    jobname = prefix

    paraset = copy.deepcopy(cmdset)
    paraset['-I'] = inputpath + samples[0] + '/' + bam
    for si in range(1,len(samples)):
        paraset['-I'] = paraset['-I'] + ' -I ' + inputpath + samples[si] + '/' + bam

    CMDs.append( cmdGenerator.formatCmd(javacmd, GATK+genotyper, paraset) )

    CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )
    jobmanager.createJob(jobname, CMDs, outpath = outputpath, outfn = jobname)
    
    return jobmanager


def picardQC(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, bam, prefix = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'bam', 'prefix'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))

    samples = cmdset.pop('sample')

    javacmd = 'java -Xmx%dg -jar'%(int(mem.replace('G',''))-2)

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    metrics = {'CollectRnaSeqMetrics.jar': 'RnaSeq', 'CollectMultipleMetrics.jar': '', 'EstimateLibraryComplexity.jar': 'Lib', 'CollectGcBiasMetrics.jar': 'GC'}
    metrickeys = ['CollectRnaSeqMetrics.jar', 'CollectMultipleMetrics.jar', 'EstimateLibraryComplexity.jar', 'CollectGcBiasMetrics.jar']
    for sample in samples:
        jobname = prefix + '_' + sample
        allcmds = []
        allcmds.append(cmdGenerator.formatCmd('source ~/libraries/setup_seqtools'))
    
        paraset = copy.deepcopy(cmdset)
        if bam == '=sample':
            paraset['INPUT'] = '=%s'%(inputpath+sample)
        else:
            paraset['INPUT'] = '=%s/%s'%(inputpath+sample, bam)
        paraset['TMP_DIR'] = paraset['TMP_DIR'] + prefix + '_' + sample + '/'
        cmdGenerator.checkPath(paraset['TMP_DIR'].strip('='), create=createpath)

        for metric in metrickeys:
            if 'MultipleMetrics' in metric:
                paraset['OUTPUT'] = '=%s'%(outputpath + sample + metrics[metric])
            else:
                paraset['OUTPUT'] = '=%s.txt'%(outputpath + sample + '.' + metrics[metric])
            paraset['CHART_OUTPUT'] = '%s'%(paraset['OUTPUT'].replace('.txt', '.pdf'))
            paraset['SUMMARY_OUTPUT'] = '%s'%(paraset['OUTPUT'].replace('.txt', '.summary.txt'))

            #filter out parameters that are not supported
            metricparaset = configRobot.validParas(paraset, availParas[metric])
            allcmds.append(cmdGenerator.formatCmd(javacmd, programpath + metric, metricparaset))

        allcmds.append(cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)))
        allcmds.append(cmdGenerator.formatCmd('rm -Rf', paraset['TMP_DIR'].strip('=')))
        jobmanager.createJob(jobname, allcmds, outpath = outputpath, outfn = jobname)
    return jobmanager

    
    
def RNASeQC(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time = configRobot.popParas(cmdset, ['cmd', 'mem', 'time'])
    samples, bam, prefix = configRobot.popParas(cmdset, ['sample', 'bam', 'prefix'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))
    if '-Djava.io.tmpdir' in cmdset.keys():
        javacmd = 'java ' + '-Djava.io.tmpdir' + cmdGenerator.checkPath(cmdset.pop('-Djava.io.tmpdir'))
    else:
        javacmd = 'java'
    javacmd = javacmd + ' -Xmx%dg -jar'%(int(mem.replace('G',''))-2)


    #need to generate -s, -o
    #generate a temperory file
    samplestr = '"Sample ID\\tBam File\\tNotes\\n'
    for sample in samples:
        samplestr = samplestr + '%s\\t%s\\t%s\\n'%(sample, inputpath + sample + '/' + bam, sample)
    samplestr = samplestr + '"'
    samplefile = '%s.samples'%prefix

    cmdset['-s'] = samplefile
    cmdset['-o'] = inputpath + 'RNA-SeQC/'

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    setupcmd = cmdGenerator.formatCmd('source ~/libraries/setup_seqtools')
    createsamplefile = cmdGenerator.formatCmd('echo', samplestr, '>', samplefile)
    removesamplefile = cmdGenerator.formatCmd('rm -f', samplefile)
    mvscriptcmd = cmdGenerator.formatCmd('mv %s%s %s'%(prefix, jobmanager.ext, cmdset['-o']))
    qccmd = cmdGenerator.formatCmd(javacmd, programpath+cmd, cmdset)

    cmdGenerator.checkPath(cmdset['-o'], create=createpath)
    jobmanager.createJob(prefix, [setupcmd, createsamplefile, qccmd, removesamplefile, mvscriptcmd], outfn = prefix, outpath = cmdset['-o'], trackcmd=False)
    return jobmanager


def RSeQC(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time = configRobot.popParas(cmdset, ['cmd', 'mem', 'time'])
    samples, bam, prefix = configRobot.popParas(cmdset, ['sample', 'bam', 'prefix'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    programs = ['inner_distance.py', 'junction_annotation.py', 'junction_saturation.py', 'read_GC.py', 'read_duplication.py']

    for sample in samples:
        jobname = prefix + '_' + sample
        CMDs = []
        CMDs.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )
        for prog in programs:
            paraset = copy.deepcopy(cmdset)
            paraset['-i'] = inputpath + sample + '/' + bam
            paraset['-o'] = outputpath + sample + '.%s'%(prog.replace('.py', ''))
            paraset = configRobot.validParas(paraset, availParas[prog])
            if '-o' not in paraset.keys():
                paraset['>'] = outputpath + sample + '.%s'%(prog.replace('.py', ''))                            

            CMDs.append( cmdGenerator.formatCmd('python', programpath+prog, paraset) )
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )
        jobmanager.createJob(jobname, CMDs, outpath = outputpath, outfn = jobname)
    return jobmanager
            

###main
def main(cmd = '', module = '', config = 'projconfig.txt', runmode = 'test'):
    global availParas
    cmdsets = configRobot.readConfig(config)
    setnames = cmdsets.keys()
    availParas = getAvailableParas()

    for setname in setnames:
        if setname != module and cmdsets[setname]['cmd'] != cmd:
            continue

        jobmanager = ''


        if 'RNA-SeQC' in cmdsets[setname]['cmd']:
            jobmanager = RNASeQC(cmdsets[setname], runmode)
        else:
            jobmanager = eval( '%s(cmdsets[setname], runmode)'%cmdsets[setname]['cmd'] )

        if runmode == 'run' and jobmanager != '':
            submitted, skipped = jobmanager.submitJob()
            jobmanager.removeJobFn(status='skipped')
            del jobmanager
            jobmanager = ''


if __name__ == '__main__':
    main( **dict( map(lambda(l): l.split('='), argv[1:]) ) )

    #usage:
    #python rnaseq-pipeline.py <option>=<val>
    #options:
    #  cmd = cmd defined in each set in configure file
    #  module = setname defined in configure file
    #  runmode = {'run', 'test'}; 'test' will only generate script files; 'run' will generate script files and submit the jobs
    #
    #  cmd and module are used to specify which command and/or which module(set) to run
    #





