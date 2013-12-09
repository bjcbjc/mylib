import argparse
#import traceback
from sys import exit
from os import system, popen

def parseLogFile(logfn, maxlinenum=100):
    #line by line, in case logfn is large
    hostname = ''
    tmpdir = []
    with open(logfn, 'r') as log:
        n = 0
        for line in log:
            n = n + 1
            if n > maxlinenum: 
                break
            if 'hostname:' in line:
                hostname = line.strip().replace('hostname:', '')
            if 'TMPJOBDIR=' in line:
                tmpdir.append( line.strip().replace('TMPJOBDIR=',''))
    return hostname, tmpdir
    
# set up parameters for the script
argp = argparse.ArgumentParser(prog='moveFileAfterJob.py')
argp.add_argument('-cp', help='copy instead of move files', action='store_true')
argp.add_argument('-fromdir', type=str, metavar='directory', default='',help='directory to move from')
argp.add_argument('-host', type=str, default='', metavar='hostname', help='hostname of the directory')
argp.add_argument('-pattern', type=str, metavar='pattern', default ='*', help='file or subdirectory pattern for moving/copying')
argp.add_argument('-log', type=str, metavar='job_log_files', help='job log files to search for directories to move/copy', default='')
argp.add_argument('-silent', help='do not print messages; default is in non-silent mode', action='store_true')
argp.add_argument('-todir', type=str, required=True, metavar='directory', default='',help='directory to move to')


args = argp.parse_known_args()
if isinstance(args, tuple):
    passon = ' '.join(args[1])
    args = args[0]

dirfile = ''
if args.fromdir != '':
    dirfile += args.fromdir
    if dirfile[-1] != '/': dirfile += '/'
if args.pattern != '':
    if args.pattern != '*' or dirfile == '':
        dirfile += args.pattern
if dirfile == '/':
    print 'I do not think you want to move/copy everything under /'
    exit()


sshcmd = 'ssh {hostname} '
if args.cp:
    mvcmd = 'cp {dirfile} {todir}'
    action = 'copy'
else:
    mvcmd = 'mv {dirfile} {todir}'
    action = 'move'

if args.log != '':
    if '*' or '?' in args.log:
        f = popen('ls %s'%args.log)
        logfns = f.read().split()
        f.close()
    else:
        logfns = [args.log]
    for jfn in logfns:
        hostname, tmpdir = parseLogFile(jfn)
        
        if not args.silent:
            print '%s files created by job %s'%(action,jfn)
            print '\tHost=%s, Dir to %s from=%s, target Dir=%s'%(hostname, action, tmpdir, args.todir)
            if dirfile != '*' and dirfile != '':
                print '\tFile pattern=%s'%(dirfile)
            else:
                dirfile = ''

        for d in tmpdir:
            if d[-1] != '/': d += '/'
            cmd = (sshcmd + mvcmd).format(hostname=hostname, dirfile=d+dirfile, todir=args.todir)
            if not args.silent:
                print '\tCmd: %s ... '%cmd,
            exitcode = system(cmd)
            if not args.silent:
                print 'exitcode=%s'%(exitcode)
else:
    if dirfile == '*' or 'scratch' not in dirfile:
        print 'I don''t think you want to %s %s'%(action,dirfile)
        exit()

    cmd = ''
    if args.host != '':
        cmd += sshcmd.format(hostname=args.host)
    cmd += mvcmd.format(dirfile=dirfile, todir=args.todir)
    
    hostname = args.host
    if hostname == '': hostname = 'localhost'

    if not args.silent:
        print '%s files'%(action)
        print '\tHost=%s, Dir to %s from=%s, target Dir=%s'%(hostname, action, dirfile, args.todir)
        print '\tCmd: %s ...'%(cmd),
    exitcode = system(cmd)
    if not args.silent:
        print 'exitcode=%s'%(exitcode)
    




    
