

#classes to manipulate jobs on cluster

from os import system
from sys import exit
import os.path
import re

class job:
    def __init__(self, fn, cmd, outfn='', outpath='', sgeopt=[], mem='4G',time='8::', status='',sgeJob=True):
        #mem and time are required; therefore separate from sgeopt
        self.fn = fn
        self.cmd = cmd
        self.outfn = outfn
        self.outpath = outpath
        self.sgeopt = sgeopt
        self.mem = mem
        self.time = time
        self.status = status #scripted, submitted, skipped
        self.fnExist = False
        self.outfnExist = False
        self.sgeJob = sgeJob
        self.checkExists()
        return

    def checkExists(self):
        if os.path.exists(self.fn):
            self.fnExist = True
        if os.path.exists(self.outpath + self.outfn):
            self.outfnExist = True
        return

class jobManager:
    def __init__(self, ext='.job', outext='.log', shell='bash', mem='4G', time='8::', overwrite=False):
        self.ext = ext
        self.outext = outext
        self.shell = shell
        self.jobs = []
        self.defaultmem = mem
        self.defaulttime = time
        self.overwrite = overwrite
        if self.overwrite == 'False': 
            self.overwrite = False
        if self.overwrite == 'True':
            self.overwrite = True
        return

    def createJob(self, fn, cmd, outfn='', outpath='/nethome/bjchen/', sgeopt=[], mem='', time='', trackcmd=True, tracklen=100, quiet=False, sgeJob=True, runOnServer='', toShell=[], runThru=False):
        #check parameters
        if type(cmd) is not list: cmd = [cmd]
        if outpath[-1] != '/': outpath = outpath + '/'
        if fn[-4:] != self.ext: fn = fn + self.ext
        if len(sgeopt) > 0 and type(sgeopt) != type([]): sgeopt = [sgeopt]
        if outfn != '' and self.outext not in outfn: outfn = outfn + self.outext
        if outfn == '':
            outpathfn = outpath + fn.replace(self.ext, '') + self.outext
        else:
            outpathfn = outpath + outfn
        if mem == '': mem = self.defaultmem
        if time == '': time = self.defaulttime
        if type(trackcmd) is str:
            trackcmd = trackcmd == 'True'
        if type(runThru) is str:
            runThru = runThru == 'True'


        if self.checkExists(fn):
            if self.overwrite:
                if not quiet:
                    print 'JobFactory.createJob: Job %s exists, OVERWRITE'%fn
            else:
                if not quiet:
                    print 'JobFactory.createJob: Job %s exists, SKIP'%fn
                return
        if self.checkExists(outpathfn):
            if self.overwrite:
                if not quiet:
                    print 'JobFactory.createJob: Log %s exists, OVERWRITE'%outpathfn
            else:
                if not quiet:
                    print 'JobFactory.createJob: Log %s exists, SKIP'%outpathfn
                return

        #write script file
        f = open( fn, 'w')
        f.write('#!/bin/%s\n'%(self.shell))
        #f.write('#$ -S /bin/%s\n'%(self.shell))

        cmdtxt = ''
        if sgeJob:
            f.write('#$ -cwd\n')
            f.write('#$ -j y\n')
            #f.write('#$ -m ea\n')
            #f.write('#$ -M bjchen@nygenome.org\n')
            #f.write('#$ -l h_vmem=%s,h_rt=%s\n'%(mem, time))
            f.write('#$ -l h_vmem=%s\n'%(mem))
            #f.write('#$ -l mem=%s,time=%s\n'%(mem, time))
            #f.write('#$ -l mem_free=%s'%(mem))
            for i in sgeopt:
                f.write('#$ %s\n'%i)

            f.write('#$ -o %s\n'%outpathfn)
        else:  
            f.write('#stdout and stderr is directed to %s'%outpathfn)
            cmdtxt = 'exec &>%s\n'%outpathfn

        f.write('\n')
        if not runThru:
            f.write('set -e\n')
        f.write('echo hostname: `hostname`\n')
        if sgeJob:
            f.write('echo jobID: $JOB_ID\n')
        if len(toShell) > 0:
            for line in toShell:
                f.write(line + '\n')
        
        cmdlineidx = []
        noncmdkey = ['echo', 'printf', 'export', 'if', 'fi', 'set ', 'shopt','alias']
        for lineidx in range( len(cmd) ):
            if any( [ cmd[lineidx].cmdstr.find(x) == 0 for x in noncmdkey ] ) or cmd[lineidx].cmdstr == '':
                continue
            cmdlineidx.append(lineidx)

        totalcmd = len(cmdlineidx)
        for lineidx in range(len(cmd)):
            cmdstr = cmd[lineidx].cmdstr
            logstr = cmd[lineidx].logstr
            if lineidx in cmdlineidx:
                if logstr != '':
                    cmdtxt = cmdtxt + ('printf "\\n======= jobCMD %d/%d: %s\\n\\n"\n'%(cmdlineidx.index(lineidx)+1, totalcmd, logstr.replace('"','""').replace('\1','\\1').replace('\2','\\2'))).replace('%', '%%')
            cmdtxt = cmdtxt + '%s\n\n'%cmdstr
            if trackcmd and lineidx in cmdlineidx:
                if logstr != '':
                    cmdtxt = cmdtxt + ('printf "finished %s\n\n"\n\n'%logstr).replace('%','%%')#[:min(len(line),tracklen)]
        cmdtxt = cmdtxt + 'echo Finish %s\n'%fn

        if runOnServer == '':
            f.write(cmdtxt)
        else:
            cmdtxt = re.sub('\n+','\n', cmdtxt)
            f.write('hostname=`hostname`\n\n')
            f.write('if [ $hostname = "%s" ]; then\n'%(runOnServer))
            f.write('\t%s'%(cmdtxt.replace('\n','\n\t').strip('\t')))
            f.write('else\n')
            f.write('\tssh %s "%s"\n'%(runOnServer, cmdtxt.replace('\n',';').replace('\\','\\\\').replace('"','\\"')))
            f.write('fi\n')
        
        f.close()

        #create a job object to keep track
        j = job(fn, cmd, outfn, outpath, sgeopt, mem, time, status='scripted', sgeJob=sgeJob)
        self.jobs.append(j)
        return 

    def submitJob(self, overwrite='', fnlist=[]):
        if overwrite == '':
            overwrite = self.overwrite
        if len(fnlist) == 0: fnlist = map(lambda(l): l.fn, self.jobs)

        alljobs = map(lambda(l): l.fn, self.jobs)
        submitted = []
        skipped = []
        for fi in range(len(fnlist)):
            f = fnlist[fi]

            if f not in alljobs: continue

            job = self.jobs[ alljobs.index(f) ] #reference
            if not overwrite:
                #delete existing log file if it os.path.exists
                if os.path.exists(job.outpath + job.outfn):
                    print 'jobFactory.submitJob: JOB %s: outfile %s exists'%(f, job.outpath + job.outfn)
                    print 'jobFactory.submitJob: Skip its submission.'
                    skipped.append(f)
                    job.status = 'skipped'
                else:
                    if job.sgeJob:
                        system('qsub %s'%f)
                    else:
                        system('bash %s'%f)
                    submitted.append(f)
                    job.status = 'submitted'
            else:
                if os.path.exists(job.outpath + job.outfn):
                    print 'jobFactory.submitJob: rm -f %s'%(job.outpath+job.outfn)
                    system('rm -f %s'%(job.outpath+job.outfn))
                if job.sgeJob:
                    system('qsub %s'%f)
                else:
                    system('nohup bash %s >/dev/null 2>&1 &'%f)
                submitted.append(f)
                job.status = 'submitted'
        return submitted, skipped

    def checkExists(self, fn = ''):
        #if fn is specified, only check for the file and return if exists or not
        #otherwise, check all jobs in the object and write the obj attributes
        if fn != '':
            return os.path.exists(fn)

        for job in self.jobs:
            job.checkExists()
        return

    def removeJobFn(self, fnlist=[], status=''):
        #delete the script file and remove the job from the job array
        if len(fnlist) == 0 and status == '':
            print 'jobFactory.removeJobFn: specify either fnlist or status'
            return
        elif len(fnlist) > 0 and status != '':
            print 'jobFactory.removeJobFn: specify either fnlist or status'
            return
        alljobs = map(lambda(l): l.fn, self.jobs)
        if len(fnlist) > 0:
            jobindex = map(lambda(a): alljobs.index(a) if a in alljobs else -1, fnlist)
            jobindex.sort(reverse=True)
            for i in jobindex:
                system('rm -f %s'%(self.jobs[i].fn))
                del self.jobs[i]

        else:
            for i in range(len(self.jobs)-1, -1, -1):
                if self.jobs[i].status == status:
                    system('rm -f %s'%(self.jobs[i].fn))
                    del self.jobs[i]
        return
