

#classes to manipulate jobs on cluster

from os import system
from sys import exit
import os.path


class job:
    def __init__(self, fn, cmd, outfn='', outpath='', sgeopt=[], mem='4G',time='8::', status=''):
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
        self.checkExists()
        return

    def checkExists(self):
        if os.path.exists(self.fn):
            self.fnExist = True
        if os.path.exists(self.outpath + self.outfn):
            self.outfnExist = True
        return

class jobManager:
    def __init__(self, ext='.job', outext='.stdout', shell='tcsh', mem='4G', time='8::', overwrite=False):
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

    def createJob(self, fn, cmd, outfn='', outpath='/ifs/scratch/c2b2/dp_lab/bc2252/', sgeopt=[], mem='', time='', trackcmd=True, tracklen=100, quiet=False):
        #check parameters
        if type(cmd) != type([]): cmd = [cmd]
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
        f.write('#$ -S /bin/%s\n'%(self.shell))
        f.write('#$ -cwd\n')
        f.write('#$ -j y\n')
        f.write('#$ -l mem=%s,time=%s\n'%(mem, time))
        for i in sgeopt:
            f.write('#$ %s\n'%i)

        f.write('#$ -o %s\n'%outpathfn)
        f.write('\n')
        for line in cmd:
            f.write('%s\n\n'%line)
            if trackcmd:
                f.write('echo finished %s\n\n'%line[:min(len(line),tracklen)])
        f.write('echo Finish %s\n'%fn)
        f.close()

        #create a job object to keep track
        j = job(fn, cmd, outfn, outpath, sgeopt, mem, time, status='scripted')
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
                    system('qsub %s'%f)
                    submitted.append(f)
                    job.status = 'submitted'
            else:
                if os.path.exists(job.outpath + job.outfn):
                    print 'jobFactory.submitJob: rm -f %s'%(job.outpath+job.outfn)
                    system('rm -f %s'%(job.outpath+job.outfn))
                system('qsub %s'%f)
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
