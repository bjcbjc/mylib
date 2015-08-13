
import subprocess
import traceback
import re
import collections
import shlex
import random
import string


class SamtoolsPileup(object):
    excludePattern = re.compile('(\^.)|(\$)')
    indelOccurenceProg = re.compile('([\+-][0-9]+)')
    def __init__(self, samtools=None, stranded=False, maxDepth=10000, minBaseQ=13, minMQ=0):
        if samtools is None:
            self.samtools = '/nfs/sw/samtools/samtools-1.1/samtools '
        else:
            self.samtools = samtools + ' '
        self.stranded = stranded
        self.args = []
        self.passOnArg = ' -B -q {minMQ} -Q {minBaseQ} -d {maxDepth} '.format(
            minMQ=minMQ, minBaseQ=minBaseQ, maxDepth=maxDepth)
        self.tempFile = []
        self.pipe = []

    def __del__(self):
        self._clean()

    def _clean_pipe(self):
        for i in xrange(len(self.pipe)):
            f = self.pipe.pop(0)
            f.close()

    def _clean(self):
        self._clean_pipe()
        for i in xrange(len(self.tempFile)):
            fn = self.tempFile.pop(0)
            subprocess.call(['rm', '-f', fn])

    def chrm_format(self, bam):
        sampipe = subprocess.Popen(shlex.split(self.samtools + ' view -H ' + bam), stdout=subprocess.PIPE)
        samheader = sampipe.stdout.read().split('\n')
        sampipe.stdout.close()
        samheader = filter(lambda(l): l[:3]=='@SQ', samheader)
        chrms = map(lambda(l): l.split()[1].replace('SN:',''), samheader)
        if 'chr' in chrms[0]:
            return 'chr'
        else:
            return ''

    def pileup_by_r(self, bam, chrm, start, end, exclude=''):
        # call this when the list of queries is short (<800 makers for example)
        chrm = self._check_loc_format(chrm)
        start = self._check_loc_format(start)
        end = self._check_loc_format(end)
        baseCounterTable = collections.defaultdict(collections.Counter)
        # loop over regions
        cmd = self.samtools + ' mpileup ' + self.passOnArg + ' -r {chrm}:{start}-{end} ' + bam
        for c, s, e in zip(chrm, start, end):
            try:
                sampipe = subprocess.Popen(shlex.split(cmd.format(chrm=c, start=s, end=e)),
                                           shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                pipe = [sampipe.stdout, sampipe.stderr]
                baseCounterTable.update(self.pileup_count(pipe[0], exclude))
            except:
                self._clean()
                print c, s, e
                traceback.print_exc()
                raise
            finally:
                self._clean_pipe()
        return baseCounterTable

    def pileup_by_l(self, bam, loc, exclude=''):
        # call this when the list of queries is long (thousands for example)
        if type(loc) is str:
            locFile = loc
        else:
            locFile = self._make_temp_file(loc)
            self.tempFile.append(locFile)
        cmd = shlex.split(self.samtools + ' mpileup ' + self.passOnArg + ' -l {locFile} '.format(locFile=locFile) + bam)
        try:
            # print cmd
            sampipe = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            pipe = [sampipe.stdout, sampipe.stderr]
            self.pipe.extend(pipe)
            return self.pileup_count(pipe[0], exclude)
        except:
            traceback.print_exc()
            raise
        # finally:
        #     self._clean()

    def pileup_count(self, samout, exclude=''):
        #pileup special chars in read_bases: ^: start of read, followed by an additional char for mapping qual
        # $: end of read
        # > or < suggests non-coding region
        if len(exclude) > 0:
            excludeProg = re.compile('[%s]+'%exclude)
        else:
            excludeProg = re.compile('')
        baseCounterTable = collections.defaultdict(collections.Counter)
        for line in samout:
            line = line.split() #chr, pos, ref, #read, read_bases, read_qual
            if len(line) > 4: #if has reads
                key = '_'.join(line[:2])
                line[4] = self.excludePattern.sub('', line[4])
                line[4] = excludeProg.sub('', line[4])
                if not self.stranded:
                    line[4] = line[4].upper().replace('<','>')

                indels, line[4] = self._matchIndels(line[4])
                #count
                baseCount = collections.Counter( line[4] )   #dict with counts and letters
                baseCount.update(indels)
                baseCounterTable[key] = baseCount
        return baseCounterTable

    @staticmethod
    def _matchIndels(baseStr):
        """ return
                indels: dict( indel, count_str)
                newBaseStr: baseStr without indels
        """
        occurrences = SamtoolsPileup.indelOccurenceProg.findall(baseStr)
        indels = collections.Counter()
        newBaseStr = baseStr
        for n in set(occurrences):
            indelprog = re.compile('\\%s[ATCGNatcgn]{%s}'%(n,n.strip('+-')))
            indelStrings = indelprog.findall(baseStr)
            newBaseStr = indelprog.sub('', newBaseStr)
            indels.update( collections.Counter(indelStrings) )
        return indels, newBaseStr

    @staticmethod
    def _check_loc_format(data):
        t = type(data)
        if t is not list and t is not tuple:
            data = [data]
        if type(data[0]) is not str:
            data = map(str, data)
        return data

    @staticmethod
    def _make_temp_file(loc):
        elementType = type(loc[0])
        if elementType is not list and elementType is not tuple:
            loc = [loc]
        fn = ''.join(random.choice(string.ascii_letters + string.digits) for x in xrange(10)) + '.pileuploc.txt'
        with open(fn, 'w') as f:
            for chrm, pos in loc:
                f.write('%s\t%s\n'%(chrm, pos))
        return fn



