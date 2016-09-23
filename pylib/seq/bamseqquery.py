
from __future__ import division
import subprocess
import traceback
import collections
import shlex
import re
import tempfile

__author__ = 'bjchen'

class SamtoolsView(object):
    def __init__(self, samtools=None, stranded=False, minMQ=0, skipFlag=3852):
        if samtools is None:
            self.samtools = '/nfs/sw/samtools/samtools-1.1/samtools '
        else:
            self.samtools = samtools + ' '
        self.stranded = stranded
        self.args = []
        self.passOnArg = ' -q {minMQ} -F {skipFlag} '.format(minMQ=minMQ, skipFlag=skipFlag)
        self.tempFile = []
        self.pipe = []
        self.cigarReduceD = re.compile('[NP]')
        self.cigarReduceM = re.compile('[X\=]')
        self.cigarParser = re.compile('(\d+)([MIDS])')

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
        samPipe = subprocess.Popen(shlex.split(self.samtools + ' view -H ' + bam), stdout=subprocess.PIPE)
        samHeader = samPipe.stdout.read().split('\n')
        samPipe.stdout.close()
        samHeader = filter(lambda(l): l[:3]=='@SQ', samHeader)
        chrms = map(lambda(l): l.split()[1].replace('SN:',''), samHeader)
        if 'chr' in chrms[0]:
            return 'chr'
        else:
            return ''

    # def query_read(self, bam, chrm, start, countSeq):
    #     """
    #
    #     Wrapper function; can be used for a single position (indel) or lists of positions
    #
    #     :param bam: bam file
    #     :param chrm: single or list of chrm
    #     :param start: single or list of pos
    #     :param end: single or list of pos
    #     :param countSeq: single or list of list of sequences that we want to count
    #     :return:
    #     """
    #
    #     chrm = self._check_loc_format(chrm)
    #     start = self._check_loc_format(start, int)
    #     countSeq = self._check_seqlist_format(countSeq)
    #
    #     counterTable = collections.defaultdict(collections.Counter)
    #     # loop over regions
    #     for c, s, seqList in zip(chrm, start, countSeq):
    #         key = (c, s)
    #         counterTable[key] = self.count_seq(bam, c, s, seqList)
    #     return counterTable

    def count_seq(self, bam, chrm, start, seqList, maxMismatch, insertion, ignoreEmptyMatchedSeq=False):
        """

        :param bam: string
        :param chrm: string
        :param start: int
        :param seqList: list/tuple of string
        :param ignoreEmptyMatchedSeq: if the matched sequenced (by position) is all '-', ignore counting;
            otherwise it would be counted as incompatible. This is likely due to intronic (spliced) regions
        :return: an instance of collections.Counter: key:seq, value:#reads support
        """

        if insertion:
            treatSAs = 'I'
        else:
            treatSAs = 'D'
        queryCounter = collections.Counter() #key: seq, value: count
        cmd = self.samtools + ' view ' + self.passOnArg + bam + ' {chrm}:{start}-{end} '
        cmd = cmd.format(chrm=chrm, start=start, end=start)
        minDiffDist = self._min_diff_seq(seqList)
        if minDiffDist is None:
            print 'seqList:', seqList
            raise ValueError('cannot find minimum sequence distinguishing seq-list\n')
        try:
            samPipe = subprocess.Popen(shlex.split(cmd), shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # pipe = [samPipe.stdout, samPipe.stderr]
            pipe = samPipe.communicate()
            # print chrm, start
            # if chrm == 'chr4' and start == 167658613:
            #     print pipe[1]
            #     print 'out'
            #     print pipe[0].split('\n')
            #     print cmd
            #     print pipe[1]
            #     exit()
            if 'random alignment retrieval only works for indexed' in pipe[1]:
                raise ValueError('index of bam %s not found'%bam)

            for line in pipe[0].split('\n'):
                line = line.split('\t') #
                if len(line) > 9: #if has reads
                    read = line[9]
                    readStart = int(line[3])

                    # process cigar
                    cigarStr = self.cigarReduceD.sub('D', line[5]) #treat N as D
                    cigarStr = self.cigarReduceM.sub('M', cigarStr) # treat =/X as M
                    cigar = [ [int(n), op, ''] for n, op in self.cigarParser.findall(cigarStr) ] #list of list [num, op, readseq]

                    # S can be I or D...
                    #skip the first S; we want extended seq support
                    if cigar[0][1] == 'S':
                        read = read[cigar[0][0]:] #clip the read
                        cigar = cigar[1:] # clip cigar

                    #treat the last S as I;
                    if cigar[-1][1] == 'S':
                        cigar[-1][1] = treatSAs
                        if not insertion:
                            read = read[:-cigar[-1][0]]
                            #cigar = cigar[:-1]

                    #push seq to cigar, cigar should only have M/I/D
                    cigar = self._matchSeqWithCigar(read, cigar)

                    # if query deletion 'del_seq' at position 'start'

                    bestMatch = [-1, None, None]
                    validQuery = False
                    for query in seqList:
                        qLen = len(query)
                        matchedCigar, matchedSeq = self._query2(cigar, readStart, start, qLen)
                        matchedLen = len(matchedCigar)
                        numMismatch = self._num_mismatch(query[:min(matchedLen, qLen)],  matchedSeq)
                        if matchedLen >= minDiffDist:
                            if ignoreEmptyMatchedSeq:
                                uniqChar = set(matchedSeq)
                                if len(uniqChar) == 1 and '-' in uniqChar:
                                    validQuery = False
                                else:
                                    validQuery = True
                            else:
                                validQuery = True
                        if matchedLen >= minDiffDist and numMismatch <= min(matchedLen-3, maxMismatch):
                            score = (matchedLen - numMismatch) / matchedLen
                            if score > bestMatch[0]:
                                bestMatch = [score, query, matchedSeq, numMismatch, matchedLen]
                        # if query[:min(matchedLen, qLen)] == matchedSeq and matchedLen >= minDiffDist:
                        # if self._num_mismatch(query[:min(matchedLen, qLen)],  matchedSeq) <= maxMismatch and matchedLen >= minDiffDist:
                        #     queryCounter[query] += 1
                        #     print query, matchedSeq, self._num_mismatch(query[:min(matchedLen, qLen)],  matchedSeq), matchedLen, line
                        # print query, matchedLen, numMismatch, matchedCigar, matchedSeq, bestMatch
                    if bestMatch[1] is not None:
                        queryCounter[bestMatch[1]] += 1
                        # print bestMatch, line
                    elif validQuery: #read not compatible and the extracted sequence for compared is not ambiguous
                        queryCounter['incompatible'] += 1
        except:
            self._clean()
            print chrm, start
            traceback.print_exc()
            raise
        finally:
            self._clean_pipe()

        return queryCounter


    @staticmethod
    def _num_mismatch(seq1, seq2):
        return sum([x != y for x, y in zip(seq1, seq2)])

    @staticmethod
    def _min_diff_seq_v0(seqList):
        """

        :param seqList: list of sequence
        :return: minimum length of sequence (from 0, first character) that differ all sequences in the list
        """
        mlen = min(map(len, seqList))
        mdist = None
        for i in xrange(mlen):
            seq = [ x[i] for x in seqList]
            if len(seq) == len(set(seq)): #all unique
                mdist = i + 1
                break
        return mdist

    @staticmethod
    def _min_diff_seq(seqList):
        """

        :param seqList: list of sequence
        :return: minimum length of sequence (from 0, first character) that differ all sequences in the list
        """
        mlen = min(map(len, seqList))
        mdist = mlen
        for i in xrange(mlen):
            seq = [ x[:i+1] for x in seqList]
            if len(seq) == len(set(seq)): #all unique
                mdist = i + 1
                break
        return mdist

    @staticmethod
    def _check_seqlist_format(data):
        """
        make sure data are nested list; eg. [['AAA','AAC'],['AAA']]
        :param data:
        :return:
        """
        t = type(data)
        if t is not list and t is not tuple:
            print 'incorrect seq list format'
            raise TypeError()
        if type(data[0]) is str:
            data = [data]
        return data

    @staticmethod
    def _check_loc_format(data, dataType=str):
        t = type(data)
        if t is not list and t is not tuple:
            data = [data]
        if type(data[0]) is not dataType:
            data = map(dataType, data)
        return data

    @staticmethod
    def _cumsum(data, offset):
        newdata = data[:]
        newdata[0] += offset
        for i in xrange(1, len(newdata)):
            newdata[i] = newdata[i-1] + newdata[i]
        return newdata

    @staticmethod
    def _matchSeqWithCigar(read, cigar):
        # break reads according to cigar
        for idx in xrange(len(cigar)): #op should have only MID
            if cigar[idx][1] in 'MI':
                cigar[idx][2] = read[:cigar[idx][0]]
                read = read[cigar[idx][0]:]
            elif cigar[idx][1] == 'D':
                cigar[idx][2] = ''
            else:
                raise ValueError('op %s in cigar, only expect MID'%cigar[idx][1])
        if len(read) > 0:
            raise ValueError('not all seq assigned to cigar: %d left; '%len(read), cigar)
        return cigar

    @staticmethod
    def _query(cigar, pos, length):
        """
        position and length imply the genomic location (minus read start)

        :param cigar: should have only M/I/D
        :param pos: position is the genomic position relative to read's first base (the first base as 1)
        :param length: length after pos, including pos
        :return:
        """

        #first find cigar op
        #need cumsum and loop from the first anyway
        endPoint = pos + length - 1
        matchedCigar = []
        matchedSeq = []
        currentEnd = 0
        lengthToFill = length
        for n, op, seq in cigar:
            currentEnd += n
            if pos <= currentEnd:
                if endPoint <= currentEnd:
                    matchedCigar.append(op * lengthToFill)
                    if op == 'D':
                        matchedSeq.append('-' * lengthToFill)
                    else:
                        matchedSeq.append(seq[pos - currentEnd + n -1 : endPoint - currentEnd + n])
                    break
                else:
                    matchedCigar.append(op * (currentEnd - pos + 1))
                    if op == 'D':
                        matchedSeq.append('-' * (currentEnd - pos + 1))
                    else:
                        matchedSeq.append(seq[pos - currentEnd + n -1 : ])
                    lengthToFill = lengthToFill - currentEnd + pos - 1
                    pos = currentEnd + 1

        matchedCigar = ''.join(matchedCigar)
        matchedSeq = ''.join(matchedSeq)
        matchedLen = len(matchedSeq)
        if len(matchedCigar)  != matchedLen:
            raise ValueError('len(matchedCigar)=%d, len(matchedSeq)=%d'%(len(matchedCigar), matchedLen))
        if matchedLen < length:
            # this should only happen if requested length is longer than cigar/read data
            if endPoint - currentEnd != length - matchedLen:
                raise  ValueError('len(matchedCigar)=%d, length=%d, endPoint=%d, currentEnd=%d'%
                                  (matchedLen, length, endPoint, currentEnd))
        elif matchedLen > length:
            raise ValueError('len(matchedCigar)=%d, length=%d'%(matchedLen, length))


        return matchedCigar, matchedSeq


    @staticmethod
    def _query2(cigar, readStart, pos, length):
        """
        position and length imply the genomic location (minus read start)

        :param cigar: should have only M/I/D
        :param readStart: genomic position of read start
        :param pos: position is the genomic position
        :param length: length after pos, including pos
        :return:
        """

        endPoint = pos + length - 1
        matchedCigar = []
        matchedSeq = []
        currentEnd = readStart - 1
        lengthToFill = length
        startFound = False

        #linear scan; shouldn't be too many cigar segments anyway
        for n, op, seq in cigar:
            if startFound: #count I, because we are asking for 'length' of seq
                # after we find the first segment containing the region of interests, we count all segments because we
                # are interested in getting a segment of sequence out from that position
                currentEnd += n
            elif op != 'I': # don't count I for genomic location
                # currentEnd is the end genomic position of the current segment of cigar
                currentEnd += n
            if pos <= currentEnd:
                if endPoint <= currentEnd:
                    matchedCigar.append(op * lengthToFill)
                    if op == 'D':
                        matchedSeq.append('-' * lengthToFill)
                    else:
                        matchedSeq.append(seq[pos - currentEnd + n -1 : endPoint - currentEnd + n])
                    break
                else:
                    matchedCigar.append(op * (currentEnd - pos + 1))
                    if op == 'D':
                        matchedSeq.append('-' * (currentEnd - pos + 1))
                    else:
                        matchedSeq.append(seq[pos - currentEnd + n -1 : ])
                    lengthToFill = lengthToFill - currentEnd + pos - 1
                    pos = currentEnd + 1
                startFound = True

        matchedCigar = ''.join(matchedCigar)
        matchedSeq = ''.join(matchedSeq)

        # error check
        matchedLen = len(matchedSeq)
        if len(matchedCigar)  != matchedLen:
            raise ValueError('len(matchedCigar)=%d, len(matchedSeq)=%d'%(len(matchedCigar), matchedLen))
        if matchedLen < length:
            # this should only happen if requested length is longer than cigar/read data
            if endPoint - currentEnd != length - matchedLen:
                raise  ValueError('len(matchedCigar)=%d, length=%d, endPoint=%d, currentEnd=%d'%
                                  (matchedLen, length, endPoint, currentEnd))
        elif matchedLen > length:
            raise ValueError('len(matchedCigar)=%d, length=%d'%(matchedLen, length))


        return matchedCigar, matchedSeq
