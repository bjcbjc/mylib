
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
    curset = {}
    #multivalpara = {}
    #multispecification = [] 
    for l in t:
        if l[0] in setnames: #new set
            if len(curset) > 0:
                settb[curset_name] = curset
                #multivalpara[curset_name] = multispecification
            curset = {}
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
                paraval = map(lambda(a): a.strip(), paraval.split(':'))
                paraval = map(int, paraval)
                paraval = map(lambda(a): '%s'%a, range(paraval[0], paraval[2]+1, paraval[1]))
                
            if paraname in curset.keys():
                if type(curset[paraname]) != type([]):
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
