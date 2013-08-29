
from os import popen, system
from sys import exit
import os.path
import vcf


def readCoverageStats():
    tb = {}
    t = open('/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/QC/bwa/COV/realign/coverage.stats').readlines()
    t.pop(0)
    t = map(lambda(l):l.strip().split(), t)
    for s, m, sd in t:
        tb[s] = [float(m), float(sd)]
    return tb

def checkduptask(scrfn):
    if len(scrfn) > 0:
        print 'multiple tasks'
        exit()

def packscr(scrfn):
    fn = []
    mem = []
    time = []
    for s in scrfn:
        if type(s) == type([]):
            fn.append(s[0])
            if len(s) >= 2: mem.append(s[1])
            if len(s) == 3: time.append(s[2])
        else:
            fn.append(s)
            mem.append('')
            time.append('')
    return fn, mem, time

def submit(fns, M, T, mem, time):
    for fi in range(len(fns)):
        f = fns[fi]
        if M[fi] == '': m = mem
        else: m = M[fi]
        if T[fi] == '': t = time
        else: t = T[fi]
        #print 'qsub -l mem=%s,time=%s %s'%(m, t, f)
        #delete existing log file if it os.path.exists
        txt = open(f).readlines()
        logfn = filter(lambda(l): '#$ -o' in l, txt)
        if len(logfn) == 1:
            logfn = logfn[0].strip('\n').strip('#$ -o ')
            if os.path.exists(logfn):
                system('rm -f %s'%logfn)
        system('qsub -l mem=%s,time=%s %s'%(m, t, f))

def genRunFile(fn, cmd, outfn='', logpath='/ifs/scratch/c2b2/dp_lab/bc2252/'):
    if logpath[-1] != '/': logpath = logpath + '/'

    f = open( fn + '.scr', 'w')
    f.write('#!/bin/tcsh\n')
    f.write('#$ -S /bin/tcsh\n')
    f.write('#$ -cwd\n')
    f.write('#$ -j y\n')
    
    if outfn == '': #output name corresponds to fn
        f.write('#$ -o %s%s.out\n'%(logpath, fn))
    else: #output to the specified name
        f.write('#$ -o %s%s.out\n'%(logpath, outfn))
    f.write('\n')

    f.write('source setup_seqtools\n')
    f.write('%s\n'%cmd)
    f.write('echo Finish %s'%fn)
    f.close()
    return fn + '.scr' #return script file name

def getAllSamBam(ext, path='./', exclude=''):
    f = popen('ls %s*%s'%(path,ext))
    t = f.read().split()
    f.close()
    if type(exclude) == type('string'):
        exclude = [exclude]
    for ex in exclude:
        if ex != '':
            t = filter(lambda(a): ex not in a, t)
    t = filter(lambda(a): a.split('/')[-1][0] != '9', t)
    return t

def getSampleName(fn):
    return fn.split('/')[-1].split('.')[0]


def removeDotExt(flist, ext):
    flist = map(lambda(l):l.replace('%s'%ext, ''), flist)
    return flist

def sam2bam(ext):
    flist = getAllSamBam(ext)
    flist = removeDotExt(flist, ext)
    for f in flist:
        system('qsub -l mem=2G,time=4:: sam2bam.scr %s'%f)

def bamindex(ext, path):
    flist = getAllSamBam(ext, path)
    for f in flist:
        system('qsub -l mem=2G,time=2:: bamindex.scr %s'%f)

def bwaalign(lanes):
    straintb = {'973': 'S288C', '916': 'Sigma'}
    for lane in lanes:
        f = popen('ls %s/*.FD.*.txt'%(lane))
        samples = f.read().split()
        f.close()
        for sample in samples:
            for s in straintb.keys():
                if s in sample:
                    samplename = sample.split('.')[2]
                    
                    cmd = 'bwa aln reference/%sseq/%s.fasta %s > BWA/%s.FD.sai\n'%(straintb[s], straintb[s], sample, samplename)
                    cmd = cmd + 'bwa aln reference/%sseq/%s.fasta %s > BWA/%s.RV.sai\n'%(straintb[s], straintb[s], sample.replace('FD', 'RV'), samplename)
                    cmd = cmd + 'bwa sampe -r \'@RG\\tID:%s\\tSM:%s\\tLB:lib\\tPL:illumina\' -n 1 reference/%sseq/%s.fasta BWA/%s.FD.sai BWA/%s.RV.sai %s %s | samtools view -bS - | samtools sort - BWA/%s\n'%(samplename, samplename, straintb[s], straintb[s], samplename, samplename, sample, sample.replace('FD','RV'), samplename)
                    scrfn.append( genRunFile(samplename, cmd) )
    return scrfn
        


def prepareGATK(ext, bampath, excludebam):
    flist = getAllSamBam(ext, bampath, excludebam)
    flist = removeDotExt(flist, ext)

    markdup = 'java -Xmx4g -jar ~/SeqTool/picard/MarkDuplicates.jar TMP_DIR=/ifs/scratch/c2b2/dp_lab/bc2252/TMP/ INPUT=%s OUTPUT=%s METRICS_FILE=%s.metrics ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT \n'
    indexbam = 'samtools index %s \n'
    createTg = 'java -Djava.io.tmpdir=/ifs/scratch/c2b2/dp_lab/bc2252/TMP/ -Xmx5g -jar ~/SeqTool/GATK/GenomeAnalysisTK.jar -I %s -R reference/%sseq/%s.fasta -T RealignerTargetCreator -o %s.intervals \n'
    align = 'java -Djava.io.tmpdir=/ifs/scratch/c2b2/dp_lab/bc2252/TMP/ -Xmx5g -jar ~/SeqTool/GATK/GenomeAnalysisTK.jar -I %s -R reference/%sseq/%s.fasta -T IndelRealigner -targetIntervals %s.intervals  -o %s \n'
    cleanup = 'rm -f %s.intervals \n'

    f = popen('ls %s*.intervals'%(bampath))
    t = removeDotExt(f.read().split(), '.intervals')
    f.close()
    existed = []
    for a in t: existed.append(getSampleName(a))

    scrfn = []
    
    for f in flist:
        if '973' in f: ref = 'S288C'
        else: ref = 'Sigma'

        mdupbam = f + '.mdup' + ext
        sample = getSampleName(f)
        
        if sample in existed:
            continue

        cmd = ''
        cmd = markdup%(f + ext, mdupbam, f + '.mdup')
        cmd = cmd + indexbam%(mdupbam)
        cmd = cmd + createTg%(mdupbam, ref, ref, f)
        cmd = cmd + align%(mdupbam, ref, ref, f, f + '.mdup.realign' + ext)
        cmd = cmd + indexbam%(f + '.mdup.realign' + ext)
        #cmd = cmd + cleanup%(f)
        scrfn.append( genRunFile(sample, cmd) )

    return scrfn
    
def rawVar_GATK(ext, bampath, excludebam, logpath, outpath, aggregate=False):

    def checkcomplete(logpath, sample):
        checkkey = ['Actual calls made', 'Total runtime', 'filtered', 'BadMateFilter', 'DuplicateReadFilter', 'UnmappedReadFilter']
        #Check if previous run os.path.exists. If the log file os.path.exists, check if the previous run has completed. If not completed, generate the script
        previousrun = False
        if os.path.exists(logpath + sample + '.out'):
            previousrun = True
            log = open(logpath + sample + '.out').readlines()
            complete = True
            for i in range(len(checkkey)):
                if checkkey[i] not in log[0-len(checkkey)+i]:
                    complete = False
                    break
        else:
            complete = False
        return complete, previousrun


    flist = getAllSamBam(ext, bampath, excludebam)
    flist = removeDotExt(flist, ext)

    scrfn = []

    javacmd = 'java -Djava.io.tmpdir=/ifs/scratch/c2b2/dp_lab/bc2252/TMP/ -Xmx%dg -jar '
    genotyper = '~/SeqTool/GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -glm BOTH -stand_call_conf 30 -stand_emit_conf 30 -R reference/%sseq/%s.fasta -o %s.vcf '

    #multiple sample call
    #groupsamples = {'S288C':'', 'Sigma':''}
    groupdrug = ['4NQO','caffeine','rapa']
    groupstrain = ['916-5','916-6','973-5','973-6']

    if not aggregate:
        for f in flist:
            sample = getSampleName(f)

            if '973' in sample: 
                ref = 'S288C'
                #groupsamples[ref] = groupsamples[ref] + '-I %s '%(f + ext)
            else: 
                ref = 'Sigma'
                #groupsamples[ref] = groupsamples[ref] + '-I %s '%(f + ext)
                
            complete, previousrun = checkcomplete(logpath, sample)
            if not complete:
                if previousrun:
                    system('rm -f %s.vcf %s.vcf.idx'%(outpath+sample, outpath+sample))
                cmd = javacmd%(4) + genotyper%(ref, ref, outpath + sample) + '-I %s \n'%(f + ext)
                scrfn.append( genRunFile(sample, cmd) )
    else:
        for s in groupstrain:
            if '973' in s:
                ref = 'S288C'
                anc = filter(lambda(l): 'ancestor' in l and '973' in l, flist)
            else:
                ref = 'Sigma'
                anc = filter(lambda(l): 'ancestor' in l and '916' in l, flist)
            for d in groupdrug:
                samparg = '-I %s '%(anc[0] + ext)
                line = s + '-' + d
                samp = filter(lambda(l): line in l, flist)
                for f in samp:
                    samparg = samparg + '-I %s '%(f + ext)
                complete, previousrun = checkcomplete(logpath, 'S'+line)
                if not complete:
                    if previousrun:
                        system('rm -f %s.vcf %s.vcf.idx'%(outpath + line, outpath + line))
                    scrfn.append( [genRunFile( 'S%s'%line, javacmd%(10) + genotyper%(ref, ref, outpath + line) + samparg + '\n', '', logpath), '12G', '120::'] )

    
    return scrfn

def rawVar_SAM(ext, bampath, excludebam, logpath, outpath, aggregate=False):

    flist = getAllSamBam(ext, bampath, excludebam)
    flist = removeDotExt(flist, ext)

    scrfn = []
    groupdrug = ['4NQO','caffeine','rapa']
    groupstrain = ['916-5','916-6','973-5','973-6']

    samcall = 'samtools mpileup -uf reference/%sseq/%s.fasta %s | bcftools view -vcg - > %s.vcf \n'

    if not aggregate:
        for f in flist:
            sample = getSampleName(f)
            if sample != 'FP-973-6-caffeine-B': continue
            if '973' in sample: 
                ref = 'S288C'
                cmd = samcall%(ref, ref, f + ext, outpath + sample)
            else: #a bug in samtools mpileup seems to complain about chr1.1 format; need to change reheader the bam and change the format in reference genome as well
                ref = 'Sigma'
                cmd = 'samtools view -H %s > %s.header\n'%(f + ext, f)
                cmd = cmd + 'python changeSigmaChrFormat.py %s.header dot2dash\n'%(f)
                cmd = cmd + 'samtools reheader %s.header %s > %s.reheader.bam \n'%(f, f + ext, f)
                cmd = cmd + samcall%(ref, ref + '.formpileup', f + '.reheader.bam', outpath + sample)
                cmd = cmd + 'python changeSigmaChrFormat.py %s.vcf dash2dot\n'%(outpath + sample)
                cmd = cmd + 'rm -f %s.header %s.reheader.bam\n'%(f, f)
            
            scrfn.append( genRunFile(sample, cmd, '', logpath) )
    else:
        for s in groupstrain:
            if '973' in s:
                ref = 'S288C'
                anc = filter(lambda(l): 'ancestor' in l and '973' in l, flist)
                for d in groupdrug:
                    samparg = '%s '%(anc[0] + ext)
                    line = s + '-' + d
                    samp = filter(lambda(l): line in l, flist)
                    for f in samp:
                        samparg = samparg + '%s '%(f + ext)
                    cmd = samcall%(ref, ref, samparg, outpath + line)            
                    scrfn.append( [genRunFile( 'S%s'%line, cmd, '', logpath), '12G', '120::'] )
            else:
                ref = 'Sigma'
                anc = filter(lambda(l): 'ancestor' in l and '916' in l, flist)
                for d in groupdrug:
                    samplist = [anc[0]]
                    line = s + '-' + d
                    samplist.extend(filter(lambda(l): line in l, flist))
                    samparg = ''
                    cmd = ''
                    for f in samplist:
                        cmd = cmd + 'samtools view -H %s > %s.header\n'%(f + ext, f)
                        cmd = cmd + 'python changeSigmaChrFormat.py %s.header dot2dash\n'%(f)
                        cmd = cmd + 'samtools reheader %s.header %s > %s.reheader.%s.bam \n'%(f, f + ext, f, line)
                        samparg = samparg + '%s.reheader.%s.bam '%(f, line)
                    
                    cmd = cmd + samcall%(ref, ref + '.formpileup', samparg, outpath + line)
                    cmd = cmd + 'python changeSigmaChrFormat.py %s.vcf dash2dot\n'%(outpath + line)
                    cmd = cmd + 'rm -f %s.header %s\n'%(f, samparg)
                                
                    scrfn.append( [genRunFile( 'S%s'%line, cmd, '', logpath), '12G', '120::'] )

    return scrfn

def rawVar_Freebayes(ext, bampath, excludebam, logpath, outpath, aggregate=False):

    flist = getAllSamBam(ext, bampath, excludebam)
    flist = removeDotExt(flist, ext)

    groupdrug = ['4NQO','caffeine','rapa']
    groupstrain = ['916-5','916-6','973-5','973-6']

    scrfn = []
    samcall = 'freebayes -p 1 -f reference/%sseq/%s.fasta -v %s.vcf '

    if not aggregate:
        for f in flist:
            sample = getSampleName(f)
            if sample != 'FP-973-6-caffeine-B': continue
            if '973' in sample: 
                ref = 'S288C'
            else: 
                ref = 'Sigma'

            cmd = samcall%(ref, ref, outpath + sample) + '-b %s'%(f + ext)
            scrfn.append( genRunFile(sample, cmd, '', logpath) )
    else:
        for s in groupstrain:
            if '973' in s:
                ref = 'S288C'
                anc = filter(lambda(l): 'ancestor' in l and '973' in l, flist)
            else:
                ref = 'Sigma'
                anc = filter(lambda(l): 'ancestor' in l and '916' in l, flist)
            for d in groupdrug:
                samparg = '-b %s '%(anc[0] + ext)
                line = s + '-' + d
                samp = filter(lambda(l): line in l, flist)
                for f in samp:
                    samparg = samparg + '-b %s '%(f + ext)
                cmd = samcall%(ref, ref, outpath + line) + samparg
                scrfn.append( [genRunFile( 'S%s'%line, cmd, '', logpath), '8G', '120::'] )

    return scrfn

def diffAncestor(ext, excludevcf, vcfpath, logpath, outpath, foutkey):
    def checkcomplete(logpath, sample):
        complete = False
        if os.path.exists(logpath + sample + '.out'):
            t = open(logpath + sample + '.out').readlines()
            if 'Total runtime' not in t[-1]:
                complete = False
            else:
                complete = True
        return complete

    flist = getAllSamBam(ext, vcfpath, excludevcf)
    flist = removeDotExt(flist, ext)

    scrfn = []

    javacmd = 'java -Djava.io.tmpdir=/ifs/scratch/c2b2/dp_lab/bc2252/TMP/ -Xmx3g -jar ~/SeqTool/GATK/GenomeAnalysisTK.jar -R reference/%sseq/%s.fasta -T VariantFiltration -o %s.%s.vcf --variant %s --mask %s --maskName Ancestor --filterExpression %s --filterName highDP '

    covtb = readCoverageStats()

    for f in flist:
        if 'ancestor' in f: continue
        sample = getSampleName(f)

        if '973' in sample: 
            ref = 'S288C'
            mask = 'BWA/ancestor973.cmb.vcf'
        else: 
            ref = 'Sigma'
            mask = 'BWA/ancestor916.cmb.vcf'

        dp_thres = covtb[sample][0] + 5.0 * covtb[sample][1]
        dp_filter = '"DP > %.1f" '%(dp_thres)

        if not checkcomplete(logpath, sample):
            cmd = javacmd%(ref, ref, outpath + sample, foutkey, f + ext, mask, dp_filter)
            scrfn.append( genRunFile(sample, cmd, '', logpath) )

    return scrfn


def combineVariants(ext, excludevcf, vcfpaths, logpath, outpath):
    def checkcomplete(logpath, sample):
        complete = False
        if os.path.exists(logpath + sample + '.out'):
            t = open(logpath + sample + '.out').readlines()
            if 'Total runtime' not in t[-1]:
                complete = False
            else:
                complete = True
        return complete
    
    flist = getAllSamBam(ext, vcfpaths[0], excludevcf)
    flist = removeDotExt(flist, ext)

    scrfn = []

    javacmd = 'java -Djava.io.tmpdir=/ifs/scratch/c2b2/dp_lab/bc2252/TMP/ -Xmx3g -jar ~/SeqTool/GATK/GenomeAnalysisTK.jar -R reference/%sseq/%s.fasta -T CombineVariants -o %s -genotypeMergeOptions UNIQUIFY -setKey set '


    for i in range(len(vcfpaths)):
        if vcfpaths[i][-1] != '/':
            vcfpaths[i] = vcfpaths[i] + '/'
    varname = map(lambda(l):l.strip('/').split('/')[-1], vcfpaths)

    for f in flist:
        sample = getSampleName(f)

        if '973' in sample:
            ref = 'S288C'
        else:
            ref = 'Sigma'

        if not checkcomplete(logpath, sample):
            varstr = ''
            for i in range(len(vcfpaths)):
                varstr = varstr + '--variant:%s %s '%(varname[i], vcfpaths[i] + sample + ext)

            cmd = javacmd%(ref, ref, outpath + sample + '.cmb.vcf') + varstr
            scrfn.append( genRunFile(sample, cmd, '', logpath) )
    return scrfn



def replaceQUAL(ext, excludevcf, vcfpath, sourcevcfpaths, sourceext, logpath, outpath, outext):
    #GATK's combineVariants mix up QUAL.. this function takes QUAL from priority order in sourcevcfpaths
    def replaceFunc(fn, sourcevcfpaths, sourceext, outpath, outext):
        for i in range(len(sourcevcfpaths)):
            if sourcevcfpaths[i][-1] != '/': sourcevcfpaths[i] = sourcevcfpaths[i] + '/'
        sourcename = map(lambda(l): l.strip('/').split('/')[-1], sourcevcfpaths)
        
        sourcedata = {}
        sample = getSampleName(fn)
        for i in range(len(sourcevcfpaths)):
            vcffile = vcf.VCFReader(open(sourcevcfpaths[i] + sample + sourceext))
            sourcedata[sourcename[i]] = []
            for r in vcffile:
                sourcedata[sourcename[i]].append( [r.CHROM, '%s'%r.POS, r.ALT, '%s'%r.QUAL] )
        txt = open(fn).readlines()
        header = filter(lambda(l):l[0] == '#', txt)
        txt = filter(lambda(l):l[0] != '#', txt)

        for li in range(len(txt)):
            record = txt[li].split('\t')
            found = False
            repstr = []
            for s in sourcename:
                match = filter(lambda(l): l[0]==record[0] and l[1]==record[1] and len(set(tuple(l[2])).intersection(set(record[4].split(',')))) > 0, sourcedata[s])
                if len(match) > 0:
                    found = True
                    A = ''
                    for fieldidx in range(5):
                        A =  A + record[fieldidx] + '\t'
                    
                    repstr.append([A+record[5], A+match[0][3]])
                    repstr.append([record[7], record[7] + ';qualSource=%s'%s])
                    break
            if not found:
                print 'Error, not found record in', fn, 'record', record
                exit()
            for A, B in repstr:
                txt[li] = txt[li].replace(A, B)
                    
        addHeader = False
        fout = open(outpath + sample + outext, 'w')
        for l in header:
            if '##INFO=<ID=' in l and not addHeader:
                fout.write('##INFO=<ID=qualSource,Number=1,Type=String,Description="QUAL Source">\n')
                addHeader = True
            fout.write('%s'%l)
        for l in txt:
            fout.write('%s'%l)
        fout.close()

    flist = getAllSamBam(ext, vcfpath, excludevcf)
    flist = removeDotExt(flist, ext)
    
    for f in flist:
        sample = getSampleName(f)
    
        replaceFunc(f + ext, sourcevcfpaths, sourceext, outpath, outext)
    return


def addSnpEff(ext, excludevcf, vcfpath, logpath, outpath, snpeffstatpath):
    flist = getAllSamBam(ext, vcfpath, excludevcf)
    flist = removeDotExt(flist, ext)

    scrfn = []

    javacmd = 'java -Djava.io.tmpdir=/ifs/scratch/c2b2/dp_lab/bc2252/TMP/ -Xmx1g -jar ~/SeqTool/snpEff/snpEff.jar eff -noLog -v -ud 1000 -c /ifs/home/c2b2/dp_lab/bc2252/SeqTool/snpEff/snpEff.config -i vcf -o vcf -s %s.html %s %s > %s '

    repcmd = 'python genenamelookup.py replace -outfile %s %s \n\n'

    for f in flist:
        sample = getSampleName(f)
        if '973' in sample:
            ref = 'S288C'
        else:
            ref = 'Sigma'

        addkey = ext.replace('.stripsnpeff.vcf','')
        cmd = javacmd%(snpeffstatpath + sample + addkey, ref, f + ext, outpath + sample + addkey + '.snpeff.vcf')
        cmd = cmd + repcmd%(outpath + sample + addkey + '.snpeff.gene.vcf', outpath + sample + addkey + '.snpeff.vcf')
        scrfn.append( genRunFile(sample, cmd, '', logpath) )

    return scrfn


def addGATKAnnotation(ext, excludevcf, vcfpath, bampath, bamext, outpath, outext, logpath):
    def checkcomplete(logpath, sample):
        complete = False
        if os.path.exists(logpath + sample + '.out'):
            t = open(logpath + sample + '.out').readlines()
            if 'Total runtime' not in t[-4]:
                complete = False
            else:
                complete = True
        return complete

    flist = getAllSamBam(ext, vcfpath, excludevcf)
    flist = removeDotExt(flist, ext)

    scrfn = []

    javacmd = 'java -Djava.io.tmpdir=/ifs/scratch/c2b2/dp_lab/bc2252/TMP/ -Xmx3g -jar ~/SeqTool/GATK/GenomeAnalysisTK.jar -T VariantAnnotator -A MappingQualityZeroFraction -A GCContent -G StandardAnnotation -R reference/%sseq/%s.fasta -V %s -L %s -I %s -o %s \n\n'
    fixtypecmd = 'python fixGATKTypeError.py %s \n'

    for f in flist:
        sample = getSampleName(f)
    
        if '973' in sample:
            ref = 'S288C'
        else:
            ref = 'Sigma'

        if not checkcomplete(logpath, sample):
            vcffile = f + ext
            cmd = javacmd%(ref, ref, vcffile, vcffile, bampath + sample + bamext, outpath + sample + outext)
            cmd = cmd + fixtypecmd%(outpath + sample + outext)
            scrfn.append( genRunFile(sample, cmd, '', logpath) )

    return scrfn


def selectVariant(ext, excludevcf, vcfpath, logpath, outpath):
    flist = getAllSamBam(ext, vcfpath, excludevcf)
    flist = removeDotExt(flist, ext)

    scrfn = []

    #need to edit the filters in vcftoolbox.py
    pycmd = 'python vcftoolbox.py selectvariant -outfile %s %s '

    for f in flist:
        sample = getSampleName(f)
        
        cmd = pycmd%(outpath + sample + '.sel.vcf', f + ext)
        scrfn.append( genRunFile(sample, cmd, 'pyselvar', logpath) )

    return scrfn


def sortVariantAndRepName(ext, excludevcf, vcfpath, logpath, outpath):
    flist = getAllSamBam(ext, vcfpath, excludevcf)
    flist = removeDotExt(flist, ext)

    scrfn = []

    #call genenamelookup.py replace
#    repcmd = 'python genenamelookup.py replace -outfile %s %s \n'
    #call vcftoolbox.py sortByIsec 
    sortcmd = 'python vcftoolbox.py sortByIsec -outfile %s %s \n'
    #remove the tmp file
    rmcmd = 'rm -f %s \n'

    for f in flist:
        sample = getSampleName(f)
        
        cmd = ''
#        cmd = repcmd%(outpath + sample + '.rep.vcf', f + ext)
#        cmd = cmd + sortcmd%(outpath + sample + '.srtg.vcf', outpath + sample + '.rep.vcf')
        cmd = cmd + sortcmd%(outpath + sample + '.srtg.vcf', f + ext)
#        cmd = cmd + rmcmd%(outpath + sample + '.rep.vcf')
        scrfn.append( genRunFile(sample, cmd, 'pysortvar', logpath) )

    return scrfn

def selectAndSortAndRepName(ext, excludevcf, vcfpath, logpath, outpath):
    flist = getAllSamBam(ext, vcfpath, excludevcf)
    flist = removeDotExt(flist, ext)

    scrfn = []

    #need to edit the filters in vcftoolbox.py
    selcmd = 'python vcftoolbox.py selectvariant -outfile %s %s \n\n'
    #call genenamelookup.py replace
#    repcmd = 'python genenamelookup.py replace -outfile %s %s \n\n'
    #call vcftoolbox.py sortByIsec 
    sortcmd = 'python vcftoolbox.py sortByIsec -outfile %s %s \n\n'
    #remove the tmp file
    rmcmd = 'rm -f %s \n'

    for f in flist:
        sample = getSampleName(f)
        
        cmd = selcmd%(outpath + sample + '.tmp1.vcf', f + ext)
#        cmd = cmd + repcmd%(outpath + sample + '.tmp2.vcf', outpath + sample + '.tmp1.vcf')
        cmd = cmd + sortcmd%(outpath + sample + '.sel.vcf', outpath + sample + '.tmp1.vcf')
        cmd = cmd + rmcmd%(outpath + sample + '.tmp1.vcf')
        scrfn.append( genRunFile(sample, cmd, 'pysortvar', logpath) )

    return scrfn

def markVariantLocInAncestor(ext, excludevcf, vcfpath, ancpaths, logpath, outpath, outext):
    flist = getAllSamBam(ext, vcfpath, excludevcf)
    flist = removeDotExt(flist, ext)
    
    scrfn = []
    cmploccmd = 'python vcftoolbox.py compareLoc -outfile %s '

    for f in flist:
        sample = getSampleName(f)
        if '973' in sample:
            anc = 'ancestor973.vcf'
        else:
            anc = 'ancestor916.vcf'

        cmpstr = ''
        for ancpath in ancpaths:
            if ancpath[-1] != '/': ancpath = ancpath + '/'
            cmpname = ancpath.strip('/').split('/')[-1]
            cmpstr = cmpstr + ' -cmpfn %s:%s '%(cmpname, ancpath + anc)

        cmd = cmploccmd%(outpath + sample + outext) + cmpstr + ' %s'%(f + ext)
        scrfn.append( genRunFile(sample, cmd, 'pymarkancloc', logpath) )
    return scrfn


def selectVariantsProcess(ext, excludevcf, vcfpath, ancpaths, logpath, outpath, filterstr, outext):
    flist = getAllSamBam(ext, vcfpath, excludevcf)
    flist = removeDotExt(flist, ext)

    scrfn = []

    #need to edit the filters in vcftoolbox.py
    selcmd = 'python vcftoolbox.py %s -outfile %s %s \n\n'
    #call genenamelookup.py replace
#    repcmd = 'python genenamelookup.py replace -outfile %s %s \n\n'
    #call vcftoolbox.py sortByIsec 
    #sortcmd = 'python vcftoolbox.py sortByIsec -outfile %s %s \n\n'
    #call vcftoolbox.py compareLoc to mark variants
    cmploccmd = 'python vcftoolbox.py compareLoc -outfile %s '
    #remove the tmp file
    rmcmd = 'rm -f %s \n'

    for f in flist:
        sample = getSampleName(f)

        if '973' in sample:
            anc = 'ancestor973.vcf'
        else:
            anc = 'ancestor916.vcf'

        if os.path.exists(outpath + sample + outext): continue

        cmd = selcmd%(filterstr, outpath + sample + '.tmp1.vcf', f + ext)
 #       cmd = cmd + repcmd%(outpath + sample + '.tmp2.vcf', outpath + sample + '.tmp1.vcf')
        #cmd = cmd + sortcmd%(outpath + sample + '.sel.vcf', outpath + sample + '.tmp2.vcf')

        cmpstr = ''
        for ancpath in ancpaths:
            if ancpath[-1] != '/': ancpath = ancpath + '/'
            cmpname = ancpath.strip('/').split('/')[-1]
            cmpstr = cmpstr + ' -cmpfn %s:%s '%(cmpname, ancpath + anc)

        cmd = cmd + cmploccmd%(outpath + sample + outext) + cmpstr + ' %s\n\n'%(outpath + sample + '.tmp2.vcf')

        cmd = cmd + rmcmd%(outpath + sample + '.tmp*.vcf')
        scrfn.append( genRunFile(sample, cmd, 'pyselvar', logpath) )
    return scrfn





#### setting ####
bamext = '.mdup.realign.bam'
excludebam = ['.mdup.bam']
vcfext = '.flt3.stripsnpeff.vcf'
excludevcf = ['.diffanc.vcf', '.deanc.vcf','.qual.vcf', '.cmb.vcf', '.snpeff.vcf']
bampath = './BWA/'
covext = '.cov1'
ifsubmit = False
mem = '2G'
time = '1::'
QCpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/QC/bwa/picard/'
FQoutpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/QC/bwa/FQ/'
covoutpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/QC/bwa/COV/realign/'
logpath = '/ifs/scratch/c2b2/dp_lab/bc2252/'
GATK_logpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/log/bwa_GATKvar/'
GATK_outpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/BWA/GATK/'
GATK_diffanclogpath  = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/log/bwa_GATK_diffanc/'
sam_outpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/BWA/SAM/'
sam_logpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/log/bwa_SAMvar/'
sam_diffanclogpath  = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/log/bwa_SAM_diffanc/'
freebayes_outpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/BWA/FB/'
freebayes_logpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/log/bwa_FBvar/'
freebayes_diffanclogpath  = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/log/bwa_FB_diffanc/'
comblogpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/log/bwa_combine/'
comboutpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/BWA/COMB/'
snpefflogpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/log/bwa_comb_snpeff/'
snpeffstatpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/BWA/COMB/snpEff_stats/'
addAntlogpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/log/bwa_addGATKAnt/'
selpath = '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/BWA/SEL/'
#################

scrfn = []


#preapre for GATK calling
#checkduptask(scrfn)
#scrfn = prepareGATK(bamext, bampath, excludebam)
 
#GATK raw variant call
#checkduptask(scrfn)
#scrfn = rawVar_GATK(bamext, bampath, excludebam, GATK_logpath, GATK_outpath, aggregate=False) 

#SAM raw variant call
#checkduptask(scrfn)
#scrfn = rawVar_SAM(bamext, bampath, excludebam, sam_logpath, sam_outpath, aggregate=False)


#diff ancestor's variants
#checkduptask(scrfn)
#scrfn = diffAncestor(vcfext, excludevcf, GATK_outpath, GATK_diffanclogpath, GATK_outpath, 'deanc')

#checkduptask(scrfn)
#scrfn = diffAncestor(vcfext, excludevcf, sam_outpath, sam_diffanclogpath, sam_outpath, 'deanc')


#combine variants from three callers
#checkduptask(scrfn)
#scrfn = combineVariants(vcfext, excludevcf, [GATK_outpath, sam_outpath, freebayes_outpath], comblogpath, comboutpath)


#replace QUAL with priority order from GATK, SAM, FB
#replaceQUAL(vcfext, excludevcf, comboutpath, [GATK_outpath, sam_outpath, freebayes_outpath], '.deanc.vcf', '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/', comboutpath, '.qual.vcf') #single thread


#add GATK annotation(stats) to the variants
#checkduptask(scrfn)
#scrfn = addGATKAnnotation(vcfext, excludevcf, comboutpath, bampath, bamext, comboutpath, '.ant.vcf', addAntlogpath)

#select, (sort), add gene name, mark same-loc-in-ancestor variants
#checkduptask(scrfn)
#scrfn = selectVariantsProcess(vcfext, excludevcf, comboutpath, [GATK_outpath, sam_outpath, freebayes_outpath], '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/', selpath, 'selectvariant', '.lmk.vcf')

#checkduptask(scrfn)
#scrfn = selectVariantsProcess(vcfext, excludevcf, comboutpath, [GATK_outpath, sam_outpath, freebayes_outpath], '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/', selpath, 'selectvariantFLT2', '.flt2.vcf')

#checkduptask(scrfn)
#scrfn = selectVariantsProcess(vcfext, excludevcf, comboutpath, [GATK_outpath, sam_outpath, freebayes_outpath], '/ifs/scratch/c2b2/dp_lab/bc2252/Gal_Seq/', selpath, 'selectvariantGSFLT3', '.flt3.vcf')


#add snpeff info to the variants
checkduptask(scrfn)
scrfn = addSnpEff(vcfext, excludevcf, selpath, snpefflogpath, selpath, selpath+'snpEff_stats/')


scrfn, M, T = packscr(scrfn)
if ifsubmit and len(scrfn) > 0:
    submit(scrfn, M, T, mem, time)


