


import os.path
from os import system

class CMD:
    def __init__(self, cmdstr, logstr=''):
        self.cmdstr = cmdstr
        self.logstr = logstr

def checkPath(p, create=False):
    if p[-1] != '/': p = p + '/'
    if create:
        if not os.path.exists(p):
            system('mkdir %s'%p.strip('='))
    return p

def checkPathOnNode(p, tmpdir=True):    
    if tmpdir:
        cmdstr = 'if [ ! -d %s ]; then \nmkdir -p %s\necho "TMPJOBDIR=%s"\nfi'%(p,p,p)
    else:
        cmdstr = 'if [ ! -d %s ]; then \nmkdir -p %s\nfi'%(p,p)
    cmd = CMD(cmdstr)
    return cmd

def passToScript(cmdset):
    #cmdset is dictionary
    #look for 'toShell' key
    key = 'toShell'
    if key in cmdset.keys():
        return formatCmd( cmdset[key] )
    else:
        return formatCmd('')

def formatCmd(cmd, *args, **kwargs ):
    #args contains tuples of all parameters for a command
    #if args[i] is any other data types (number, string, list), it's printed directly
    #if args[i] is a dict, it's printed as 'key value'
    #order of args should follow the order of the command 'cmd'
    #
    # kwargs: logstr='', string for logging command(s); for tracking purpose
    #

    if 'logstr' in kwargs:
        logstr = kwargs['logstr']
    else:
        logstr = ''

    pipeout = ''
    if type(cmd) == type([]):
        cmdstr = '\n'.join(cmd)
    else:
        cmdstr = '%s'%cmd
    for i in args:
        if type(i) == type({}):
            if '>' in i.keys():
                pipeout = i.pop('>')
            for k in i.keys():
                if type(i[k]) == type([]) or type(i[k]) == type(()):
                    if i[k][0][0] == '=':
                        for a in i[k]:
                            cmdstr = cmdstr + ' %s=%s'%(k, a.strip('='))
                    else:
                        for a in i[k]:
                            cmdstr = cmdstr + ' %s %s'%(k, a)
                else:
                    if len(i[k]) == 0:
                        cmdstr = cmdstr + ' %s'%k
                    elif i[k][0] == '=':
                        cmdstr = cmdstr + ' %s%s'%(k, i[k])
                    else:
                        cmdstr = cmdstr + ' %s %s'%(k, i[k])
        elif type(i) == type([]) or type(i) == type(()):
            for k in i:
                cmdstr = cmdstr + ' %s'%k
        else:
            cmdstr = cmdstr + ' %s'%i
    if pipeout != '':
        cmdstr = cmdstr + ' > %s'%pipeout

    if logstr == '':
        logstr = cmdstr
    cmdobj = CMD(cmdstr, logstr)
    return cmdobj
    
    



