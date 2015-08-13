#!~/bin/python2.7

from __future__ import division
from sys import argv
from os import system, popen
import os.path
import copy
import cmdGenerator
import jobFactory
import configRobot
import itertools
import re
import string
import random
import math



def randstr(size=10, chars=string.ascii_letters + string.digits):
    return ''.join(random.choice(chars) for x in range(size))

def readSampleList(fn):
    sample = [line.strip().split('\t')[0] for line in open(fn)]
    return sample

def readInputList(fn):
    data = [line.strip().split('\t') for line in open(fn)]
    return data

def runJobs(cmdset, runmode='test'):
    jobmanager = jobFactory.jobManager(mem=cmdset['mem'], time=cmdset['time'], overwrite=cmdset.pop('overwrite'))
    if 'sampleList' in cmdset:
        cmdset['LoopOverSample'] = readSampleList(cmdset['sampleList'])
    if 'LoopOverLine' in cmdset:
        cmdset['LoopOverLine'] = readInputList(cmdset['LoopOverLine'])
    loopVarName, loopList, cmdset = configRobot.getLoopOverList( cmdset )

    if type(cmdset['call']) is not list:
        cmdset['call'] = [cmdset['call']]

    if 'njob' in cmdset:
        njob = int(cmdset['njob'])
        chunksize = int( math.ceil(len(loopList)/njob))
    else:
        chunksize = 1
    
    #unpack call
    cmdset['call_idx'] = []
    for i in xrange(len(cmdset['call'])):
        key = 'call_cmd_%d'%i
        cmdset[key] = cmdset['call'][i]
        cmdset['call_idx'].append(key)
    del cmdset['call']

    allowance = chunksize
    CMDs = []
    sampleParser = re.compile('Sample_([\w\-]+)')
    for loopValue in loopList:
        paraset = copy.deepcopy(cmdset)
        if 'randstr' in paraset:
            paraset['randstr'] = randstr(20)
        paraset.update( zip( loopVarName, loopValue))
        if 'LoopOverLine' in paraset:
            if type(paraset['LoopOverLine']) is not list:
                paraset['LoopOverLine'] = [paraset['LoopOverLine']]
            for fdIdx in xrange(len(paraset['LoopOverLine'])):
                paraset['line%d'%fdIdx] = paraset['LoopOverLine'][fdIdx]
        if 'sample' not in paraset:
            if 'LoopOverFile' in paraset:
                if 'Sample_' in paraset['LoopOverFile']:
                    paraset['sample'] = sampleParser.findall(paraset['LoopOverFile'])[0]
        if 'logpath' not in paraset: paraset['logpath'] = paraset['outputpath']
        paraset = configRobot.evalRegExp( paraset )
        paraset = configRobot.translateAllValues( paraset )
        
        if 'tmppath' in paraset:
            CMDs.append( cmdGenerator.checkPathOnNode( paraset['tmppath'] ) )
        CMDs.append( cmdGenerator.checkPathOnNode( paraset['outputpath'] ) )
        CMDs.append( cmdGenerator.checkPathOnNode( paraset['logpath'] ) )             
        for callKey in cmdset['call_idx']:
            CMDs.append( cmdGenerator.formatCmd( paraset[callKey] ) )
        allowance -= 1
        if allowance == 0:
            CMDs.append( cmdGenerator.formatCmd( 'mv {prefix}.job {logpath}'.format(**paraset) ))
            CMDs.append( cmdGenerator.formatCmd( 'mv {prefix}.log {logpath}'.format(**paraset) ))
            jobmanager.createJob(paraset['prefix'], CMDs, outpath = './', outfn = paraset['prefix'], trackcmd=paraset['trackcmd'], sgeJob=True, sgeopt=paraset['sgeopt'], toShell=paraset['toShell'], runThru=paraset['runThru'])
            CMDs = []
            allowance = chunksize
    else:
        if len(CMDs) > 0:
            CMDs.append( cmdGenerator.formatCmd( 'mv {prefix}.job {logpath}'.format(**paraset) ))
            CMDs.append( cmdGenerator.formatCmd( 'mv {prefix}.log {logpath}'.format(**paraset) ))
            jobmanager.createJob(paraset['prefix'], CMDs, outpath = './', outfn = paraset['prefix'], trackcmd=paraset['trackcmd'], sgeJob=True, sgeopt=paraset['sgeopt'], toShell=paraset['toShell'], runThru=paraset['runThru'])
            CMDs = []
            allowance = chunksize

    return jobmanager



###main
def main(module = '', runmode = 'test', config=os.path.dirname(os.path.abspath(__file__))+'/job.config.txt'):
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
