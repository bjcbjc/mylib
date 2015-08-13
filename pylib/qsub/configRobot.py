
#functions to read configure files and prepare settings for running programs
#
#Configure file contains setup, paratermeters for running programs on sets of data
#Each 'set' is defined by a '<name>:'.
#Under each set, program parameters can be set directly parameter name should follow the program's usage.
#That is, if the program allows -para1 and --para2, then to set up the values for these parameters, you should sepecify '-para1' and '--para2' (note single '-' and double '-')
#For each line in a 'set', the format is <parameter> <tab> <values>; values can be empty string if the paramter is a switch (taking no values), if multiple values are provided, they should be separated by ';'. If the paramter is a range, use start:step:end to specify the range (a list of values)
#
#Required parameters for each set: 
#  'mem', 'time' since these will specify the requirement of a job
#  'outputpath': path for output files, including a file containing messages from stdout
#  'prefix': prefix is used for output files (is possible) and job script, stdout files; it might be also used to create job specific directory (<outpath>/<prefix>/), depending on how each job is created. This is useful when a program, such as tophat, create generic output files (no option to change output files), in which case, each job needs to output the file to its own directory.
#
#Configure files can have comments. Lines starting with '#' will be ignored.
#

import glob
import itertools
import operator                
import traceback
import copy
import re
def readConfig(configfn):
    t = map(lambda(l):l.strip(), open(configfn).readlines())
    t = filter(lambda(l): len(l) > 0 and l[0] != '#', t)
    t = map(lambda(l): l.split('\t'), t)
    setnames = filter(lambda(l): len(l) == 1 and ':' in l[0], t)
    setnames = map(lambda(l): l[0], setnames)

    #some program supports specify a parameter multiple times (eg. PROGRAM=A, PROGRAM=B...)
    #store these values into list, but restore them later; indicate these parameters in multispecification


    #print setnames
    settb = {}
    defaultset = {'sgeopt':'', 'toShell':[], 'trackcmd':True, 'runThru':False}
    defaultn = len(defaultset)
    curset = copy.deepcopy( defaultset)
    #multivalpara = {}
    #multispecification = [] 
    for l in t:
        if l[0] in setnames: #new set
            if len(curset) > defaultn:
                settb[curset_name] = curset
                #multivalpara[curset_name] = multispecification
            curset = copy.deepcopy(defaultset)
            #multispecification = []
            curset_name = l[0].strip(':')
        else: #parameter line
            paraname = l[0]
            if len(l) > 1:
                paraval = l[1]
            else:
                paraval = ''
            if ';' in paraval: paraval = map(lambda(a): a.strip(), paraval.split(';'))
            if paraval.count(':') == 2 and paraname != 'time':
                original_paraval = paraval
                try:
                    paraval = map(lambda(a): a.strip(), paraval.split(':'))
                    paraval = map(int, paraval)
                    paraval = map(lambda(a): '%s'%a, range(paraval[0], paraval[2]+1, paraval[1]))
                except:
                    paraval = original_paraval
                
            if paraname in curset:
                if paraname in defaultset:
                    if type(curset[paraname]) is list:
                        curset[paraname].append(paraval)
                    else:
                        curset[paraname] = paraval
                else:
                    if type(curset[paraname]) is not list:
                        curset[paraname] = [curset[paraname]]
                    curset[paraname].append(paraval)
                #if paraname not in multispecification:
                #    multispecification.append(paraname)                
            else:
                curset[paraname] = paraval
    #append last curset
    settb[curset_name] = curset
    #multivalpara[curset_name] = multispecification
    return settb #, multivalpara

def popParas(paraset, paranames):
    #pop out values specified by paranames from paraset
    #values are popped out in order specified in paranames, in tuple format
    if type(paranames) != type([]): paranames = [paranames]
    paravals = map(lambda(a): paraset.pop(a), paranames)
    if len(paravals) == 1:
        paravals = paravals[0]
    return paravals

def validParas(paraset, paralist):
    #return a new dict that only contains keys that are in paralist
    return dict( (k, v) for k,v in paraset.items() if k in paralist)

def makeParasList(paraset, paranames):
    #check each parameter and if it's not list, make it list
    if type(paranames) != type([]): paranames = [paranames]
    for a in paranames:
        if a not in paraset.keys(): continue
        if type(paraset[a]) != type([]):
            paraset[a] = [ paraset[a] ]
    return paraset

def popLoopOver(paraset):
    loopOver = {}
    for k in paraset.keys():
        if 'LoopOver' in k:
            if type(paraset[k]) != type([]):
                loopOver[k] = [ paraset[k] ]
            else:
                loopOver[k] = paraset[k]
            del paraset[k]
            #see if there are wildcards
            if any( [ w in x for x in loopOver[k] for w in ['*', '?']] ):
                loopOver[k] = reduce(operator.add, map(glob.glob, loopOver[k]))
    return loopOver, paraset

def getLoopOverList(paraset):
    loopOver, paraset = popLoopOver(paraset)
    keys = loopOver.keys()
    iterlist = list( itertools.product( *[ loopOver[k] for k in keys]) )
    #keys will be names in cmdset that start with 'LoopOver'
    #iterlist will be the list containg tuples (same number of LoopOver variables) for looping
    #the reason to convert it to list is so that it can be used within a loop
    return keys, iterlist, paraset

def translateAllValues(paraset, maxiter=10):
    n = 0
    for k in paraset.keys():
        if type(paraset[k]) is str:
            paraset[k] = paraset[k].replace('\{', '#(').replace('\}', '#)')
    while n < maxiter:
        for k in paraset.keys():
            if type(paraset[k]) is not str: continue
            if '{' in paraset[k]:
                try:
                    paraset[k] = paraset[k].format( **paraset)
                except:
                    print k
                    print paraset[k]
                    print traceback.print_exc()
                    exit()
        invalid = [ k for k, v in paraset.iteritems() if type(v) is str if '{' in v ]
        if len(invalid) == 0:
            break
        else:
            n = n + 1
    if len(invalid) > 0:
        print 'not filled parameters:'
        print invalid
        exit()
    for k in paraset.keys():
        if type(paraset[k]) is str:
            paraset[k] = paraset[k].replace('#(', '{').replace('#)', '}')
            #paraset[k] = paraset[k].replace('#(', '{{').replace('#)', '}}')
    return paraset

def evalRegExp(paraset):
    for key, value in paraset.iteritems():
        if type(value) is str:
            if 're.' in value and '(' in value:
                #translate first
                if '{' in value:
                    value = value.format( **paraset )
                paraset[key] = eval( value )

    return paraset

def evalGlob(paraset):
    for key, value in paraset.iteritems():
        if type(value) is str:
            if 'glob.' in value and '(' in value:
                if '{' in value:
                    value = value.format( **paraset )
                paraset[key] = eval( value )
    return paraset
    
