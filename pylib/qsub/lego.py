#!~/bin/python2.7


from sys import argv
from os import system, popen
import os.path
import copy
import cmdGenerator
import jobFactory
import configRobot
import itertools


matlabcmd = '/nfs/apps/matlab/current/bin/matlab -nodisplay -r "addpath(genpath(\'/ifs/home/c2b2/dp_lab/bc2252/libraries/Matlabox\')); addpath(genpath(\'./functions\')); addpath(genpath(\'./regression\')); addpath(genpath(\'./data\')); %s; quit"\n'
matlabcmd2012 = '/nfs/apps/matlab/2012a/bin/matlab -nodisplay -r "addpath(genpath(\'/ifs/home/c2b2/dp_lab/bc2252/libraries/Matlabox\')); addpath(genpath(\'./functions\')); addpath(genpath(\'./regression\')); addpath(genpath(\'./data\')); %s; quit"\n'


def testSim(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmdset = configRobot.makeParasList(cmdset, ['njob', 'dataset', 'cumulatecount', 'costCoef', 'priormethod', 'transresidual', 'btthres', 'adaptiveSigmoid', 'splitpriornorm'])

    cmd, call, mem, time, prefix, njobs = configRobot.popParas(cmdset, ['cmd', 'call', 'mem', 'time', 'prefix', 'njob'])
    datasets, algo, noiseidx, cumulatecount, costCoef, priormethod, transresidual, btthres, adaptiveSigmoid, splitpriornorm = configRobot.popParas(cmdset, ['dataset', 'algo', 'noiseidx', 'cumulatecount', 'costCoef', 'priormethod','transresidual', 'btthres', 'adaptiveSigmoid', 'splitpriornorm'])

    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    cmdGenerator.checkPath(outputpath, create=createpath)

    iterlist = list(itertools.product( datasets, cumulatecount, noiseidx, costCoef, priormethod, transresidual, btthres, adaptiveSigmoid, splitpriornorm ) )

    for dataset, cumcount, nsidx, betacost, prmethod, transres, btthreshold, adaptSig, splitNorm in iterlist:
        if len(njobs) == 1:
            numjobs = int(njobs[0])
        else:
            numjobs = int(njobs[ datasets.index(dataset) ] )
        if cumcount == 'true': cumStr = 'Cum1'
        else: cumStr = 'Cum0'
        if prmethod == 'bootstrap': prStr = 'BT'
        elif prmethod == 'bootstraplinear': prStr = 'BL'
        elif prmethod == 'bootstrapexpected': prStr = 'BE'
        elif prmethod == 'bayes': prStr = 'BY'
        else: 
            print 'unknown prior method:', prmethod
            exit(1)
        if transres == 'true':
            if prmethod != 'bootstrap': continue
            else: resStr = 'Resd1' + '_' + btthreshold
        else: 
            resStr = 'Resd0'
            if btthreshold != '0.3': continue
            else: resStr = resStr + '_0.3'
        if adaptSig == 'true': sigStr = 'ADSig1'
        else: sigStr = 'ADSig0'
        if splitNorm == 'true': splitnormStr = 'SPN1'
        else: splitnormStr = 'SPN0'
        
        if splitNorm == 'true' and prmethod != 'bootstrapexpected': continue

        resStr = resStr.replace('.','')

        Nstr = 'N%s'%(nsidx)
        #costStr = 'P' + betacost
        costStr = ''
        cumStr = ''

        for jobidx in range(1, numjobs+1):
            fnhead = prefix + dataset + cumStr + costStr + Nstr + prStr + resStr + sigStr + splitnormStr
            jobname = prefix + dataset + cumStr + costStr + Nstr + prStr + resStr + sigStr + splitnormStr + 'J%02d'%(jobidx)                         
            CMD = []

            functionCall = "tic; "

            functionCall = functionCall + call + "(%d, %d, '%s', '%s', %s, '%s', '%s', 'cumulatecount', %s, 'costCoef', %s, 'priormethod', '%s', 'transresidual', %s, 'adaptiveSigmoid', %s, 'btthres', %s, 'splitpriornorm', %s); toc;"%(jobidx, numjobs, dataset, algo, nsidx, outputpath, fnhead, cumcount, betacost, prmethod, transres, adaptSig, btthreshold, splitNorm)

            if algo == 'lasso':
                CMD.append( cmdGenerator.formatCmd( matlabcmd2012%(functionCall) ) )
            else:
                CMD.append( cmdGenerator.formatCmd( matlabcmd%(functionCall) ) )
            CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )
            jobmanager.createJob(jobname, CMD, outpath=outputpath, outfn=jobname, trackcmd=False)
    return jobmanager



def runCV(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmdset = configRobot.makeParasList(cmdset, ['njob', 'dataset', 'pheno', 'usecn', 'useDisease', 'combineSplit', 'ylogtransform'])

    cmd, call, mem, time, prefix, njobs = configRobot.popParas(cmdset, ['cmd', 'call', 'mem', 'time', 'prefix', 'njob'])
    datasets, phenos, usecn, useDisease, combineSplit = configRobot.popParas(cmdset, ['dataset', 'pheno', 'usecn', 'useDisease', 'combineSplit'])
    ylogtransform, trainmode = configRobot.popParas(cmdset, ['ylogtransform', 'trainmode'])
    if 'useallexp' in cmdset.keys():
        useallexp = configRobot.popParas(cmdset, ['useallexp'])
    else: useallexp = 'false'
    if 'samplinghead' in cmdset.keys():
        samplinghead = configRobot.popParas(cmdset, ['samplinghead'])
    else: samplinghead = 'Samplings'
    if 'runcv' in cmdset.keys():
        runcv = range(1, int(configRobot.popParas(cmdset, ['runcv']))+1)
    else: runcv = [0]

    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    cmdGenerator.checkPath(outputpath, create=createpath)


    iterlist = list(itertools.product( datasets, phenos, usecn, combineSplit, useDisease, ylogtransform, runcv ) )

    for dataset, pheno, cnSet, splitSet, diseaseSet, ylog, cvidx in iterlist:
        if pheno in ['ACT', 'GI50'] and ylog == 'true': continue
        if len(njobs) == 1:
            numjobs = int(njobs[0])
        else:
            numjobs = int(njobs[ datasets.index(dataset) ] )
        if cnSet == 'true': cnStr = 'C1'
        else: cnStr = 'C0'
        if splitSet == 'true': splitStr = 'S1'
        else: splitStr = 'S0'
        if diseaseSet == 'true': diseaseStr = 'D1'
        else: diseaseStr = 'D0'
        if ylog == 'true': logStr = 'log'
        else: logStr = ''
        if cvidx == 0: cvStr = ''
        else: cvStr = '%02d'%(cvidx)

        if pheno == '[]': pheno = ''

        for jobidx in range(1, numjobs+1):
            fnhead = prefix + cvStr + logStr + cnStr + splitStr + diseaseStr
            jobname = prefix + cvStr + logStr + dataset + pheno + cnStr + splitStr + diseaseStr + 'J%02d'%(jobidx)                         
            CMD = []

            if dataset.lower() == 'joe':
                functionCall = "addpath(\'./Gray/data\'); tic; "
            else:
                functionCall = "addpath(\'./CCLE/data\'); tic; "

            functionCall = functionCall + call + "(%d, %d, '%s', '%s', '%s', '%s', 'usecn', %s, 'combineSplit', %s, 'useDisease', %s, 'ylogtransform', %s, 'trainmode', %s, 'sampling', '%s'); toc;"%(jobidx, numjobs, outputpath, fnhead, dataset, pheno, cnSet, splitSet, diseaseSet, ylog, trainmode, samplinghead+cvStr )

            CMD.append( cmdGenerator.formatCmd( matlabcmd2012%(functionCall) ) )
            CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )
            jobmanager.createJob(jobname, CMD, outpath=outputpath, outfn=jobname, trackcmd=False)
    return jobmanager




def runBootstrapCCLE(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmdset = configRobot.makeParasList(cmdset, ['njob', 'dataset', 'pheno', 'usecn', 'useDisease', 'combineSplit', 'priormethod', 'ylogtransform', 'useallexp', 'cumulatecount', 'predList'])

    cmd, call, mem, time, prefix, njobs = configRobot.popParas(cmdset, ['cmd', 'call', 'mem', 'time', 'prefix', 'njob'])
    iter, priormethod, algo, ylogtransform = configRobot.popParas(cmdset, ['iter', 'priormethod', 'algo', 'ylogtransform'])
    datasets, phenos, usecn, combineSplit, useDisease = configRobot.popParas(cmdset, ['dataset', 'pheno', 'usecn', 'combineSplit', 'useDisease'])
    useallexp = configRobot.popParas(cmdset, ['useallexp'])
    if 'trainmode' in cmdset.keys():
        trainmode = configRobot.popParas(cmdset, ['trainmode'])
    else: trainmode = 'false'
    if 'cumulatecount' in cmdset.keys():
        cumulatecount = configRobot.popParas(cmdset, ['cumulatecount'])
    if 'runcv' in cmdset.keys():
        runcv = range(1, int(configRobot.popParas(cmdset, ['runcv']))+1)
    else: runcv = [0]
    if 'samplinghead' in cmdset.keys():
        sampling = configRobot.popParas(cmdset, ['samplinghead'])
    else: sampling = 'Samplings'
    if 'lassoparahead' in cmdset.keys():
        lassopara = configRobot.popParas(cmdset, ['lassoparahead'])
    else: lassopara = 'lassCvPara'
    if 'predList' in cmdset.keys():
        predlist = configRobot.popParas(cmdset, ['predList'])
    else:
        predlist = ['cancerGenes']

    iterlist = list(itertools.product( datasets, phenos, usecn, combineSplit, useDisease, priormethod, ylogtransform, useallexp, cumulatecount, runcv, predlist ) )

    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    cmdGenerator.checkPath(outputpath, create=createpath)

    f = popen('ls %s*.finished'%(outputpath))
    finfns = f.read().split()
    f.close()
    finfns = map(lambda(l):l.replace(outputpath, '').replace('.finished', ''), finfns)

    for dataset, pheno, cnSet, splitSet, diseaseSet, prior, ylog, useExp, cumcount, cvidx, plist in iterlist:
        if pheno in ['ACT', 'GI50'] and ylog == 'true': continue
        if len(njobs) == 1:
            numjobs = int(njobs[0])
        else:
            numjobs = int(njobs[ datasets.index(dataset) ] )
        if cnSet == 'true': cnStr = 'C1'
        else: cnStr = 'C0'
        if splitSet == 'true': splitStr = 'S1'
        else: splitStr = 'S0'
        if diseaseSet == 'true': diseaseStr = 'D1'
        else: diseaseStr = 'D0'
        if prior == 'bootstrapnorm': priorStr = 'NRM'
        else: priorStr = ''
        if ylog == 'true': logStr = 'log'
        else: logStr = ''
        if useExp == 'true': expStr = 'AE'
        else: expStr = ''
        if cumcount == 'true': cumStr = 'Cum'
        else: cumStr = ''
        if cvidx == 0: cvStr = ''
        else: cvStr = '%02d'%(cvidx)
        if plist == 'cancerGenes': plStr = ''
        else: plStr = plist

        if pheno == '[]': pheno = ''

        for jobidx in range(1, numjobs+1):
            fnhead = prefix + cvStr + plStr + priorStr + logStr + cnStr + splitStr + diseaseStr + expStr + cumStr
            jobname = prefix + cvStr + plStr + priorStr + logStr + dataset + pheno + cnStr + splitStr + diseaseStr + expStr + cumStr + 'J%02d'%(jobidx)                         
            CMD = []
            for iteridx in range(1, int(iter)+1):
                if fnhead+dataset+'_J%02dI%02d'%(jobidx,iteridx) in finfns: continue
                if dataset.lower() == 'joe':
                    functionCall = "addpath(\'./Gray/data\'); tic; "
                elif 'ccle' in dataset.lower():
                    functionCall = "addpath(\'./CCLE/data\'); tic; "
                
                functionCall = functionCall + call + "(%d, %d, %d, '%s', '%s', '%s', '%s', 'usecn', %s, 'combineSplit', %s, 'useDisease', %s, 'algo', '%s', 'priormethod', '%s', 'ylogtransform', %s, 'trainmode', %s, 'useallexp', %s, 'cumulatecount', %s, 'sampling', '%s', 'lassopara', '%s', 'predList', '%s'); toc;"%(jobidx, numjobs, iteridx, outputpath, fnhead, dataset, pheno, cnSet, splitSet, diseaseSet, algo, prior, ylog, trainmode, useExp, cumcount, sampling+cvStr, lassopara+cvStr, plist )

                if algo == 'lasso':
                    CMD.append( cmdGenerator.formatCmd( matlabcmd2012%(functionCall) ) )
                else:
                    CMD.append( cmdGenerator.formatCmd( matlabcmd%(functionCall) ) )
            if len(CMD) > 0:
                CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )
                jobmanager.createJob(jobname, CMD, outpath=outputpath, outfn=jobname, trackcmd=False)
    return jobmanager


def runBootstrap(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, call, mem, time, prefix, datasets, phenos, njobs = configRobot.popParas(cmdset, ['cmd', 'call', 'mem', 'time', 'prefix', 'dataset', 'pheno', 'njob'])
    iter, usecn, combineSplit, useDisease = configRobot.popParas(cmdset, ['iter', 'usecn', 'combineSplit', 'useDisease'])

    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    cmdGenerator.checkPath(outputpath, create=createpath)

    if type(datasets) != type([]): datasets = [datasets]
    if type(njobs) != type([]): njobs = [njobs]
    if type(usecn) != type([]): usecn = [usecn]
    if type(combineSplit) != type(combineSplit): combineSplit = [combineSplit]

    for di in range(len(datasets)):
        numjobs = int(njobs[di])
        for cnSet in usecn:
            if cnSet == 'true': cnStr = 'C1'
            else: cnStr = 'C0'
            for splitSet in combineSplit:
                if splitSet == 'true': splitStr = 'S1'
                else: splitStr = 'S0'
                for jobidx in range(1, numjobs+1):
                    CMD = []
                    fnhead = prefix + datasets[di]  + cnStr + splitStr
                    jobname = prefix + datasets[di]  + cnStr + splitStr + 'J%02d'%(jobidx)                         
                    for iteridx in range(1, int(iter)+1):                        
                        functionCall = "addpath(\'./Gray/data\'); tic; " + call + "(%d, %d, %d, '%s', '%s', 'usecn', %s, 'combineSplit', %s, 'useDisease', %s); toc;"%(jobidx, numjobs, iteridx, outputpath, fnhead, cnSet, splitSet, useDisease )

                        CMD.append( cmdGenerator.formatCmd( matlabcmd%(functionCall) ) )
                    CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )
                    jobmanager.createJob(jobname, CMD, outpath=outputpath, outfn=jobname, trackcmd=False)
    return jobmanager


def gseaTestAnt2Pheno(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, call, mem, time, prefix, datasets, phenos, features, njobs = configRobot.popParas(cmdset, ['cmd', 'call', 'mem', 'time', 'prefix', 'dataset', 'pheno', 'feature', 'njob'])

    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    cmdGenerator.checkPath(outputpath, create=createpath)

    for di in range(len(datasets)):
        numjobs = int(njobs[di])
        if datasets[di] != 'CCLE': 
            continue
        for jobidx in range(1, numjobs+1):
            jobname = prefix + datasets[di] + '%02d'%jobidx
            functionCall = call + "(%d, %d, '%s', 'dataset', '%s', 'pheno', '%s', 'feature', '%s');"%(jobidx, numjobs, outputpath + prefix, datasets[di], phenos[di], features[di])
            CMD = []
            CMD.append( cmdGenerator.formatCmd( matlabcmd%(functionCall) ) )
            CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )
            jobmanager.createJob(jobname, CMD, outpath=outputpath, outfn=jobname, trackcmd=False)
    return jobmanager



def gmfit(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, call, mem, time, prefix, Ks, datafn = configRobot.popParas(cmdset, ['cmd', 'call', 'mem', 'time', 'prefix', 'Ks', 'datafn'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    cmdGenerator.checkPath(outputpath, create=createpath)

    if type(Ks) != type([]):
        Ks = [Ks]

    for fn in datafn:
        for k in Ks:
            jobname = prefix + fn.replace('.mat', '') + 'k%03d'%(int(k))
            functionCall = call + "('%s', %s, '%s');"%(inputpath+fn, k, outputpath + jobname + '.mat')
            CMD = []
            CMD.append( cmdGenerator.formatCmd( matlabcmd%(functionCall) ) )
            CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )
            jobmanager.createJob(jobname, CMD, outpath=outputpath, outfn=jobname, trackcmd=False)
    return jobmanager



###main
def main(module = '', runmode = 'test', config='projconfig.txt'):
    cmdsets = configRobot.readConfig(config)
    setnames = cmdsets.keys()

    for setname in setnames:
        if setname != module:
            continue

        jobmanager = ''

        jobmanager = eval( '%s(cmdsets[setname], runmode)'%cmdsets[setname]['cmd'] )

        if runmode == 'run' and jobmanager != '':
            submitted, skipped = jobmanager.submitJob()
            jobmanager.removeJobFn(status='skipped')
            del jobmanager
            jobmanager = ''



if __name__ == '__main__':
    main( **dict( map(lambda(l): l.split('='), argv[1:]) ) )
