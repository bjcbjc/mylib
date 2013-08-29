

#obtain available parameters for each program and stores them in a text file so other scripts can determine which paramters are valid

from os import popen
from subprocess import Popen, PIPE, STDOUT

programpaths = open('setup_seqtools').readlines()
programpaths = filter(lambda(l): len(l)>0, map(lambda(l):l.strip(),programpaths))

program = ['picard', 'tophat', 'RNA-SeQC', 'cufflinks', 'cuffcompare', 'cuffmerge', 'cuffdiff']
programpathgroup = ['picard', 'tophat', 'RNA-SeQC', 'cufflinks', 'cufflinks', 'cufflinks', 'cufflinks']

programpaths = filter(lambda(l): sum(map(lambda(a): a.lower() in l.lower(), programpathgroup))>0, programpaths)


paras = {}
for p in program:
    ppath = filter(lambda(l): programpathgroup[program.index(p)] in l, programpaths)
    if len(ppath) > 1: 
        print 'more than one path for %s'%p
        print 'ignore %s'%p
        continue
    ppath = ppath[0]
    ppath = filter(lambda(l): '/' in l, ppath.split())
    ppath = ppath[0].strip('"').replace(':','').replace('$PATH','')
    if ppath[-1] != '/': ppath = ppath + '/'

    if p == 'picard':
        f = popen('ls %s*.jar'%ppath)
        jars = f.read().split()
        f.close()
        for jar in jars:
            if 'sam-' in jar or 'picard-' in jar:
                continue
            f = Popen('java -Xmx1g -jar %s -H'%jar, shell=True, stdout=PIPE, stderr=STDOUT).stdout
            t = f.readlines()
            f.close()
            t = filter(lambda(l): len(l)>0, map(lambda(l):l.strip('\n'), t))
            t = filter(lambda(l): l[0] != ' ', t)
            t = map(lambda(l):l.split()[0], t)
            options = t[ t.index('Options:')+1: ]
            options = map(lambda(l): l.split('=')[0], options)
            paras[jar.split('/')[-1]] = options
    elif p == 'tophat':
        f = Popen('tophat', shell=True, stdout=PIPE, stderr=STDOUT).stdout
        t = f.readlines()
        f.close()
        t = filter(lambda(l):len(l)>0, map(lambda(l): l.strip('\n'), t))
        opt = map(lambda(l):l.split()[0], filter(lambda(l):l[0] == '-', map(lambda(l):l.strip(), t[t.index('Options:'):])))
        options = []
        for o in opt:
            if '/' in o: options.extend(o.split('/'))
            else: options.append(o)
        paras[p] = options
    elif p in ['cufflinks', 'cuffcompare', 'cuffmerge', 'cuffdiff']:
        f = Popen(p, shell=True, stdout=PIPE, stderr=STDOUT).stdout
        t = f.readlines()
        f.close()
        t = filter(lambda(l):len(l)>0, map(lambda(l): l.strip('\n'), t))
        opt = map(lambda(l): l.split()[0], filter(lambda(l):l.strip()[0]=='-' and len(l.split())>1, t))
        options = []
        for o in opt:
            if '/' in o: options.extend(o.split('/'))
            else: options.append(o)
        paras[p] = options
    elif p == 'RNA-SeQC':
        f = Popen('java -Xmx1g -jar ~/SeqTool/RNA-SeQC/RNA-SeQC_v1.1.4.jar', shell=True, stdout=PIPE, stderr=STDOUT).stdout
        t = f.readlines()
        f.close()
        options = map(lambda(l):l.split()[0], filter(lambda(l):l[0] == '-', map(lambda(l): l.strip(), t)))
        paras[p] = options




f = open('AllAvailableParas.txt','w')
for p, opts in paras.iteritems():
    f.write('%s'%p)
    for o in opts: f.write('\t%s'%o)
    f.write('\n')
f.close()
