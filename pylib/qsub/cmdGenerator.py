


import os.path
from os import system

def checkPath(p, create=False):
    if p[-1] != '/': p = p + '/'
    if create:
        if not os.path.exists(p):
            system('mkdir %s'%p.strip('='))
    return p

def checkPathOnNode(p):
    cmd = 'if [ ! -d %s ]; then \nmkdir %s\necho "TMPJOBDIR=%s"\nfi'%(p,p,p)
    return cmd

def formatCmd(cmd, *args ):
    #args contains tuples of all parameters for a command
    #if args[i] is any other data types (number, string, list), it's printed directly
    #if args[i] is a dict, it's printed as 'key value'
    #order of args should follow the order of the command 'cmd'

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
    return cmdstr
    
    



