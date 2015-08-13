
import subprocess
import shlex
import fastaAccess
import glob
import random
import string
import bamseqquery

__author__ = 'bjchen'

class IndelCounter(object):
    def __init__(self, refFn, exonFn, tmpDir, bedtools, samtools=None, minMQ=0, skipFlag=3852):
        self.refFn = refFn
        self.exonFn = exonFn
        self.ref = None
        self.refChrm = None
        self.tmpDir = tmpDir
        self.bedtools = bedtools
        self.samtools = samtools
        self.samViewEngine = bamseqquery.SamtoolsView(samtools=samtools, minMQ=minMQ, skipFlag=skipFlag)
        # find fai
        fn = glob.glob('%s.fai'%self.refFn)
        if len(fn) == 1:
            self.refFai = fastaAccess.readFai(fn[0])
        else:
            raise ValueError('%d fai found.'%(len(fn)))

    def query_ref(self, chrm, start, end):
        if self.refChrm != chrm:
            self.read_ref(chrm)
        try:
            return self.ref[start-1:end]
        except:
            print 'access outside chromosome region: %s, %d, %d'%(chrm, start, end)
            raise

    def read_ref(self, chrm):
        self.ref = fastaAccess.getChrmSeq(self.refFai, self.refFn, chrm)
        self.refChrm = chrm
        if self.ref is None: #chrm not in reference
            raise ValueError('%s not in reference'%chrm)

    def intersect_exon(self, indelRegions):
        firstLine = open(self.exonFn).readline().split()[0]
        if 'chr' in firstLine:
            chrmFormat = 'chr'
        else:
            chrmFormat = ''
        tmpBed = self.tmpDir + '/' + ''.join([random.choice(string.ascii_letters) for _ in xrange(20)]) + '.bed'
        with open(tmpBed, 'w') as f:
            for idx, reg in enumerate(indelRegions):
                reg[0] = self.correct_chrm(reg[0], chrmFormat)
                reg[1:] = map(str, reg[1:])
                f.write('\t'.join(reg) + '\t%d\n'%idx) #include index as identifier

        #call bedtools
        cmd = self.bedtools + ' intersect -a %s -b %s '%(tmpBed, self.exonFn)

        out = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE).stdout
        exonic_region = []
        for line in out:
            chrm, start, end, idx = line.strip().split('\t')
            exonic_region.append((int(idx), (chrm, int(start), int(end))))
        out.close()
        subprocess.call(shlex.split('rm -f %s'%tmpBed))
        exonic_region = list(set(exonic_region))
        exonic_region.sort()
        # if len(exonic_region) != len(set([x[0] for x in exonic_region])):
        #     print exonic_region
        #     raise ValueError('index not unique; multiple exonic regions overlapping indels\n')
        return exonic_region


    def output_indel_count_result(self, result, outputFile):
        header = ['gene', 'chrm', 'start', 'end', 'ref', 'alt', 'extStart', 'extEnd', 'extRef', 'extAlt', 'extRefCount', 'extAltCount', 'extIncompatibleCount']
        with open(outputFile, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for line in result:
                f.write('\t'.join(map(str, line)) + '\n')

    def count_indel(self, bamFile, variantFile, maxMismatch):
        # get indel regions in inputFile; [ [chrm, start, end, ref, alt], [originStart, originEnd, originRef, originAlt] ]
        indels = self.process_indel(variantFile)

        chrmFormat = self.get_chrm_format(bamFile)
        #filter to get exonic indels; [ [chrm, start, end, idx] ]
        exonicIndels = self.intersect_exon([reg[:3] for reg, dummy in indels])

        indelCount = dict()
        #now, for each exonic indel, we do counting
        for idx, reg in exonicIndels:
            chrm, start, end = reg
            chrm = self.correct_chrm(chrm, chrmFormat)
            # if chrm != '3' and start != 52443607:
            #     continue

            # because regions might be cut due to exon filtering, we need to adjust ref and alt as well
            ref0, alt0 = indels[idx][0][3:5]
            ref = ref0[ start - indels[idx][0][1] : end - indels[idx][0][1] + 1].upper()
            alt = alt0[ start - indels[idx][0][1] : end - indels[idx][0][1] + 1].upper()
            if ref == alt: #after exonic filtering, the sequences are the same, this just suggest it's outside exons
                continue
            try:
                counter = self.samViewEngine.count_seq(bamFile, chrm, start, [ref, alt], maxMismatch, '-' in ref)
            except:
                print idx, reg
                print ref0, alt0
                print ref, alt
                raise

            if len(counter) > 0:
                if idx in indelCount:
                    if counter[alt] > indelCount[idx][-2]:
                        indelCount[idx] = [start, end, ref, alt, counter[ref], counter[alt], counter['incompatible']]
                else:
                    indelCount[idx] = [start, end, ref, alt, counter[ref], counter[alt], counter['incompatible']]
            # else:
            #     print 'not found:', chrm, start, end, ref, alt
            #     print counter
            #     exit()

        result = list()
        #consolidate results
        for idx in xrange(len(indels)):
            record = [indels[idx][0][-1], indels[idx][0][0]]
            record.extend(indels[idx][1])
            if idx in indelCount:
                record.extend(indelCount[idx])
            else:
                record.extend([''] * 4 + ['0'] * 3)
            result.append(record)
        return result

    def process_indel(self, variantFile):
        """
        return list of indels; for each element, it contains three lists
        [chrm, matched_start, matched_end], [ref, alt], [original_start, original_end]]
        :param variantFile:
        :return:
        """
        indels = []
        if variantFile[-4:] == '.vcf':
            print 'not supported yet'
            exit()
            # for line in open(inputFile):
            #     if line[0] != '#':
            #         line = line.strip().split('\t')
            #         chrm, start, ref, alt = line[0], line[1], line[3], line[4]
            #         if len(ref) > 1 or len(alt) > 1: #indel
            #             indels.append( [[chrm, start], [ref, alt]] )
        elif variantFile[-4:] == '.maf':
            with open(variantFile) as f:
                header = f.readline()
                while header[0] == '#':
                    header = f.readline()
                for line in f:
                    line = line.strip().split('\t')
                    gene, chrm, start, end, ref, alt = line[0], line[4], int(line[5]), int(line[6]), line[11], line[12]
                    ref2 = ref.replace('-', '')
                    alt2 = alt.replace('-', '')
                    if len(ref2) == 1 and len(alt2) == 1:
                        continue
                    elif len(ref2) > len(alt2): #deletion, extend base for ref
                        start2, end2, ref2, alt2 = self.search_del_lookup(chrm, start, end, ref2)
                    else: #insertion, alt is the inserted seq
                        start2, end2, ref2, alt2 = self.search_ins_lookup(chrm, start, alt2)
                    indels.append([[chrm, start2, end2, ref2, alt2, gene], [start, end, ref, alt]])
        else:
            raise ValueError('unknown file format ', variantFile)
        return indels

    def search_del_lookup(self, chrm, start, end, seq):
        """
        seq is the deleted sequence on chrm:start-end; extend the regions on both ends to include actual sequences
        for lookup
        :param chrm:
        :param start:
        :param end:
        :param seq:
        :return:
        """
        n = len(seq)
        start, end = start - n, end + n
        ref = self.query_ref(chrm, start, end)
        alt = ref[:n] + '-' * n + ref[2*n:]
        if len(ref) != len(alt):
            raise ValueError('extending del look-up error; len(ref)=%d, len(alt)=%d'%(len(ref), len(alt)))
        return start, end, ref, alt

    def search_ins_lookup(self, chrm, start, seq):
        """
        seq is the sequence inserted at chrm:start; extend the region up/down stream to include neighboring sequences;
        however, because seq can be 'repeated' sequences in the genome, so extension is done to include until the inserted
        seq is not found
        :param chrm:
        :param start:
        :param seq: inserted sequence
        :return:
        """
        n = len(seq)
        insertionLoc = start
        end = start
        start -= n
        up = self.query_ref(chrm, start, start + n - 1)
        while up == seq:
            start -= n
            up = self.query_ref(chrm, start, start + n - 1)

        end += n
        down = self.query_ref(chrm, end - n + 1, end)
        while down == seq:
            end += n
            down = self.query_ref(chrm, end - n + 1, end)

        # at this point, start:start+n and end-n+1:end+1 are of len(seq) but different from seq
        # this is to take care when the insertion is repeated sequence
        # however, the position is reference based; the sequence we want to compare should be
        # reference[start-1:end+len(seq)] and reference[start-1:insertionLoc-1] + seq + reference[insertion_loc-1:end]
        ref = self.query_ref(chrm, start, end + n)
        alt = self.query_ref(chrm, start, insertionLoc-1) + seq + self.query_ref(chrm, insertionLoc, end)
        if len(ref) != len(alt):
            print chrm, start, end, n, insertionLoc, len(seq), up, down
            raise ValueError('extending ins look-up error; len(ref)=%d, len(alt)=%d'%(len(ref), len(alt)))
        return start, end + n, ref, alt

    def get_chrm_format(self, bamFile):
        cmd = self.samtools + ' view -H %s'%bamFile
        out = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE).stdout
        line = out.readline()
        line = out.readline()
        out.close()
        if 'chr' in line:
            chrmFormat = 'chr'
        else:
            chrmFormat = ''
        return chrmFormat

    @staticmethod
    def correct_chrm(chrm, chrmFormat):
        if chrmFormat == 'chr' and 'chr' not in chrm:
            chrm = 'chr' + chrm
            chrm = chrm.replace('MT', 'M')
        elif chrmFormat == '' and 'chr' in chrm:
            chrm = chrm.replace('chr', '')
            if chrm == 'M':
                chrm = 'MT'
        return chrm




if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-variantFile', type=str, required=True, help='variant file, maf')
    parser.add_argument('-bamFile', type=str, required=True, help='rna bam file')
    parser.add_argument('-outputFile', type=str, required=True, help='output file')
    parser.add_argument('-reference', type=str,
                        default='/data/NYGC/Resources/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa',
                        help='reference genome')
    parser.add_argument('-exonAnnt', type=str,
                        default='/nethome/bjchen/DATA/GenomicInfo/GENCODE/gencode.v18.annotation.exon.gtf',
                        help='exon annotation for filtering for RNA')
    parser.add_argument('-tmpDir', type=str, default='./', help='directory for temporary files')
    parser.add_argument('-bedtools', type=str,
                        default='/data/NYGC/Software/bedtools/bedtools-2.17.0/bin/bedtools',
                        help='path to bedtools')
    parser.add_argument('-samtools', type=str,
                        default='/nfs/sw/samtools/samtools-1.1/samtools',
                        help='path to samtools')
    parser.add_argument('-minMQ', type=int, default=0, help='min mapping quality for a read to be counted')
    parser.add_argument('-maxMismatch', type=int, default=2, help='maximum number of mimatches allowed')
    parser.add_argument('-skipFlag', type=int, default=3852, help='sam flag to filter reads')


    args = vars(parser.parse_args())


    counter = IndelCounter(args['reference'], args['exonAnnt'], args['tmpDir'], args['bedtools'],
                           args['samtools'], args['minMQ'], args['skipFlag'])

    result = counter.count_indel(args['bamFile'], args['variantFile'], args['maxMismatch'])
    counter.output_indel_count_result(result, args['outputFile'])

