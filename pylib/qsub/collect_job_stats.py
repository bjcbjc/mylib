import subprocess
import shlex
import sys
import os.path
import getpass
import smtplib

__author__ = 'bjchen'


class QacctParser(object):
    def __init__(self, fn = None):
        self.fn = fn
        self.keep = ['jobname', 'jobnumber', 'ru_wallclock', 'maxvmem', 'qname', 'hostname', 'slots', 'cpu']
        self.jobData = dict()

    def parse(self, fn = None):
        if fn is not None:
            if self.fn is None:
                self.fn = fn
            elif self.fn != fn:
                print 'Overwrite file for parsing: %s -> %s'%(self.fn, fn)
                self.fn = fn
        if self.fn is None:
            print 'No file to parse.'
            return

        jobdata = dict()
        for line in open(self.fn):
            if line.startswith('==='):
                if len(jobdata) == len(self.keep): #valid job
                    self.jobData[int(jobdata['jobnumber'])] = jobdata
                jobdata = dict()
            else:
                line = line.strip().split()
                if line[0] in self.keep:
                    jobdata[line[0]] = ' '.join(line[1:])
        return

    def output(self, outfn, col=None):
        if col is not None:
            if type(col) is str:
                col = [col]
            invalid = set(col).difference(self.keep)
            if len(invalid) > 0:
                raise ValueError('Invalid columns: %s'%invalid)
        else:
            col = self.keep

        ids = self.jobData.keys()
        ids.sort()
        if type(outfn) is str:
            f = open(outfn, 'w')
            openFile = True
        else:
            f = outfn
            openFile = False
        f.write('\t'.join(col) + '\n')
        for jobId in ids:
            f.write('\t'.join([self.jobData[jobId][k] for k in col]) + '\n')

        if openFile:
            f.close()
        return


class JobIdCollection(object):
    """
    improve speed by query keyword (prefix) in qacct first, so only one qacct call
    """
    def __init__(self, fn, addKeys = ('ru_wallclock', 'maxvmem', 'cpu', 'failed', 'exit_status')):
        self.fn = fn
        self.addKeys = addKeys
        self.keySet = set(self.addKeys)
        self.data = list()
        self.prefix = None
        self.user = getpass.getuser()
        self.jobStats = dict()
        self.errorJobIds = list()

        self.__read_job_id()


    def __read_job_id(self):
        """
        read all data
        :return:
        """
        for line in open(self.fn):
            line = line.strip().split()  # name, job id, submission status
            self.data.append(line)
            self.jobStats[line[1]] = dict()

        self.prefix = os.path.commonprefix([x[0] for x in self.data])
        if len(self.prefix) == 0:
            raise ValueError('cannot find common prefix among job names')

    def collect_job_stats(self, outFn = None):
        if outFn is None:
            out = sys.stdout
        else:
            out = open(outFn, 'w')

        self.get_job_stats()
        for line in self.data:
            if line[2] == 'submitted':
                jobStat = self.jobStats[line[1]]
                addData = [jobStat[key] if key in jobStat else 'NA' for key in self.addKeys]
                newLine = '\t'.join(line + addData)
                out.write(newLine + '\n')

                if len(jobStat) > 0:
                    if jobStat['failed'] != '0' or jobStat['exit_status'] != '0':
                        self.errorJobIds.append(line[1])
        if outFn is not None:
            out.close()
            if len(self.errorJobIds) > 0:
                msg = []
                for jobid in self.errorJobIds:
                    jobInfo = 'jobname: {jobname}\njobnumber: {jobnumber}\nfailed: {failed}\nexit_status{exit_status}\n\n'
                    jobInfo = jobInfo.format(self.jobStats[jobid])
                    msg.append(jobInfo)

                msg = ''.join(msg)
                with open(outFn + 'errorjobs', 'w') as out:
                    out.write(msg)
                self.email_errors(msg)
        return

    def email_errors(self, errorMsg):
        server = smtplib.SMTP('localhost')
        emailAddr = self.user + '@nygenome.org'
        msg = 'From: {emailAddr}\nTo: {emailAddr}\nSubject: Errors in jobs {prefix}*\n\n'
        msg = msg.format(emailAddr=emailAddr, prefix=self.prefix)
        problem = server.sendmail(emailAddr, emailAddr, msg + errorMsg)
        print 'email problem:', problem
        server.quit()

    def get_job_stats(self):
        p = subprocess.Popen(shlex.split('/opt/sge/bin/lx-amd64/qacct -u %s -j "%s*"' % (self.user, self.prefix)),
                             stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
        jobStat = dict()
        for line in p.stdout: #jobs are separated by ===
            if line.startswith('====='): #new job
                if len(jobStat) > 0:
                    if jobStat['jobnumber'] in self.jobStats:
                        self.jobStats[jobStat['jobnumber']] = jobStat
                jobStat = dict()
            else:
                ## we don't need to store everything!
                line = line.strip().split()
                if line[0] in self.keySet or line[0] in {'jobnumber', 'jobname'}:
                    jobStat[line[0]] = ' '.join(line[1:])
        else:
            if len(jobStat) > 0:
                if jobStat['jobnumber'] in self.jobStats:
                    self.jobStats[jobStat['jobnumber']] = jobStat
        p.stdout.close()
        return


class JobIdCollectionOld(object):
    def __init__(self, fn, addKeys = ('ru_wallclock', 'maxvmem', 'cpu')):
        self.fn = fn
        self.addKeys = addKeys

    def collect_job_stats(self, outFn = None):
        if outFn is None:
            out = sys.stdout
        else:
            out = open(outFn, 'w')

        for line in open(self.fn):
            line = line.strip().split() #name, job id, submission status
            if line[2] == 'submitted':
                job_stat = self.get_job_stats(line[1])
                addData = [job_stat[key] if key in job_stat else 'NA' for key in self.addKeys]
                newLine = '\t'.join(line + addData)
                out.write(newLine + '\n')
        if outFn is not None:
            out.close()
        return

    @staticmethod
    def get_job_stats(jobId):
        p = subprocess.Popen(shlex.split('qacct -j %s' % jobId), stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
        data = [line.strip().split() for line in p.stdout]
        data = {line[0]: ' '.join(line[1:]) for line in data}
        p.stdout.close()
        return data


if __name__ == '__main__':
    line = open(sys.argv[1]).readline()
    if line.startswith('====='):
        qp = QacctParser(sys.argv[1])
        qp.parse()
        if len(sys.argv) > 2:
            qp.output(sys.argv[2])
        else:
            qp.output(sys.stdout)
    else:
        jc = JobIdCollection(sys.argv[1])
        if len(sys.argv) > 2:
            jc.collect_job_stats(sys.argv[2])
        else:
            jc.collect_job_stats()