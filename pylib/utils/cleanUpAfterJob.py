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
argp = argparse.ArgumentParser(prog='cleanUpAfterJob.py')
argp.add_argument('-dir', type=str, metavar='directory', default='',help='directory to clean up')
argp.add_argument('-host', type=str, default='', metavar='hostname', help='hostname of the directory')
argp.add_argument('-pattern', type=str, metavar='pattern', default ='*', help='file or subdirectory pattern for deleltion')
argp.add_argument('-log', type=str, metavar='job_log_files', help='job log files to search for directories to delete', default='')
argp.add_argument('-silent', help='do not print messages; default is in non-silent mode', action='store_true')
#argp.add_argument('-onlyEmptyDir', action='store_true', help='delete only empty directories')

args = argp.parse_known_args()
if isinstance(args, tuple):
    passon = ' '.join(args[1])
    args = args[0]

dirfile = ''
if args.dir != '':
    dirfile += args.dir
    if dirfile[-1] != '/': dirfile += '/'
if args.pattern != '':
    if args.pattern != '*' or dirfile == '':
        dirfile += args.pattern
if dirfile == '/':
    print 'I do not think you want to delete everything under /'
    exit()


sshcmd = 'ssh {hostname} '
rmcmd = 'rm -rf {dirfile} '

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
            print 'Clean up job %s'%(jfn)
            print '\tHost=%s, Dir to clean=%s'%(hostname, tmpdir)
            if dirfile != '*' and dirfile != '':
                print '\tFile pattern=%s'%(dirfile)
            else:
                dirfile = ''

        for d in tmpdir:
            if d[-1] != '/': d += '/'
            cmd = (sshcmd + rmcmd).format(hostname=hostname, dirfile=d+dirfile)
            if not args.silent:
                print '\tCmd: %s ... '%cmd,
            exitcode = system(cmd)
            if not args.silent:
                print 'exitcode=%s'%(exitcode)
else:
    if dirfile == '*' or 'scratch' not in dirfile:
        print 'I don''t think you want to delete %s'%dirfile
        exit()

    cmd = ''
    if args.host != '':
        cmd += sshcmd.format(hostname=args.host)
    cmd += rmcmd.format(dirfile=dirfile)
    
    hostname = args.host
    if hostname == '': hostname = 'localhost'

    if not args.silent:
        print 'Clean up'
        print '\tHost=%s, Dir to clean=%s'%(hostname, dirfile)
        print '\tCmd: %s ...'%(cmd),
    exitcode = system(cmd)
    if not args.silent:
        print 'exitcode=%s'%(exitcode)
    




    
