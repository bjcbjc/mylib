
from __future__ import division

__author__ = 'bjchen'


def correct_chrm(chrm):
    if 'chr' in chrm:
        chrm = chrm.replace('chr', '')
        if chrm == 'M':
            chrm = 'MT'
    return chrm

def get_allele(fn):
    fileType = fn[-3:]

    if fileType == 'vcf':
        query = ['#CHROM', 'POS', 'REF', 'ALT']
    elif fileType == 'maf':
        query = ['Chromosome', 'Start_Position', 'End_Position', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']
    else:
        raise ValueError('unsupported file type %s'%fileType)

    result = dict()
    with open(fn) as f:
        line = f.readline()
        if fileType == 'vcf':
            while line[:6] != '#CHROM':
                line = f.readline()
        elif fileType == 'maf':
            while line[0] == '#':
                line = f.readline()

        header = line.strip().split('\t')
        idx = [header.index(x) for x in query ]
        for line in f:
            line = line.split('\t')
            data = [line[i] for i in idx]
            data[0] = correct_chrm(data[0])
            if fileType == 'maf':
                result['_'.join(data[:3])] = [x.upper() for x in data[3:]]
            else:
                result['_'.join([data[0], data[1], data[1]])] = [x.upper() for x in data[2:]]
    return result


class MAF(object):
    def __init__(self, fn, query = None,
                 keys = ('Chromosome', 'Start_Position', 'End_Position', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2')):
        if query is None:
            self.query = ['Gene_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']
        else:
            self.query = query
        self.data = dict()
        self.keys = list()  #small enough to keep
        with open(fn) as f:
            line = f.readline()
            while line[0] == '#':
                line = f.readline()
            header = line.strip().split('\t')
            keyIdx = [header.index(x) for x in keys]
            queryIdx = [header.index(x) for x in self.query ]
            if 'Chromosome' in self.query:
                chrIdx = self.query.index('Chromosome')
            else:
                chrIdx = None
            for line in f:
                line = line.split('\t')
                key = '_'.join([line[i] for i in keyIdx])
                data = [line[i] for i in queryIdx]
                if chrIdx is not None:
                    data[chrIdx] = correct_chrm(data[chrIdx])
                self.data[key] = dict(zip(self.query, data))
                self.keys.append(key)

    def get_ref_alt(self):
        res = dict()
        if 'Tumor_Seq_Allele1' not in self.query or 'Tumor_Seq_Allele2' not in self.query:
            raise ValueError('not allele info')
        for key, record in self.data.iteritems():
            res[key] = (record['Tumor_Seq_Allele1'], record['Tumor_Seq_Allele2'])
        return res

    def get_vaf_table(self, pileup, indel, outputFn):
        """

        :param pileup: PileupTable
        :param indel: IndelCountTable
        :param outputFn:
        :return:
        """
        indelTable = indel.get_vaf_dp()
        header = ['gene', 'chr', 'start', 'ref', 'alt', 'nRef', 'nAlt', 'nOther', 'DP', 'VAF']
        with open(outputFn, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for key in self.keys: #keep the order
                line = [self.data[key][fd] for fd in ['Gene_Symbol', 'Chromosome', 'Start_Position', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']]
                if '-' in self.data[key]['Tumor_Seq_Allele1'] or '-' in '-' in self.data[key]['Tumor_Seq_Allele2']:
                    if key in indelTable:
                        line.extend([str(indelTable[key][fd]) for fd in ['refCount', 'altCount', 'otherCount', 'dp', 'vaf']])
                    else:
                        line.extend(['0']*5)
                else:
                    poskey = '_'.join([self.data[key][fd] for fd in ['Chromosome', 'Start_Position', 'End_Position']])
                    tb = pileup.get_vaf_dp_single(poskey, self.data[key]['Tumor_Seq_Allele1'], self.data[key]['Tumor_Seq_Allele2'])
                    line.extend([str(tb[fd]) for fd in ['refCount', 'altCount', 'otherCount', 'dp', 'vaf']])
                f.write('\t'.join(line) + '\n')


class rnaVariantCount(object):
    def __init__(self, fn, keys, chrm):
        with open(fn) as f:
            self.header = f.readline().strip('\n').split('\t')
            self.keyIdx = [self.header.index(x) for x in keys]
            self.data = [line.strip('\n').split('\t') for line in f]
            chrIdx = self.header.index(chrm)
            if 'chr' in self.data[0][chrIdx]:
                for i in xrange(len(self.data)):
                    self.data[i][chrIdx] = correct_chrm(self.data[i][chrIdx])

class IndelCountTable(rnaVariantCount):
    def __init__(self, fn):
        super(IndelCountTable, self).__init__(fn, ['chrm', 'start', 'end'], 'chrm')
        self.countIdx = [self.header.index(key) for key in ('extRefCount', 'extAltCount', 'extIncompatibleCount')]

    def get_vaf_dp(self):
        res = dict()
        alleleIdx = [self.header.index(x) for x in ['ref', 'alt']]
        for record in self.data:
            key = [record[i] for i in self.keyIdx]
            key.extend([record[i] for i in alleleIdx])
            refCount, altCount, otherCount = map(int, [record[i] for i in self.countIdx])
            dp = refCount + altCount
            if dp > 0:
                vaf = altCount / dp
            else:
                vaf = 0
            res['_'.join(key)] = {'refCount': refCount, 'altCount': altCount, 'otherCount': otherCount,
                                  'dp': dp, 'vaf': vaf}
        return res


class PileupTable(rnaVariantCount):
    def __init__(self, fn):
        validNt = ('A', 'T', 'C', 'G', 'N', '*', 'a', 't', 'c', 'g','n')
        super(PileupTable, self).__init__(fn, ['chr', 'pos', 'pos'], 'chr')
        self.ntIdx = [idx for idx in xrange(len(self.header)) if self.header[idx] in validNt]
        data = self.data
        self.data = dict()
        for record in data:
            key = '_'.join(record[i] for i in self.keyIdx)
            self.data[key] = record

    def get_vaf_dp_single(self, key, ref, alt):
        upperHeader = [x.upper() for x in self.header]
        refCount, altCount, otherCount = 0, 0, 0
        if key in self.data:
            record = self.data[key]
            for i in self.ntIdx:
                if upperHeader[i] == ref:
                    refCount += int(record[i])
                elif upperHeader[i] == alt:
                    altCount += int(record[i])
                else:
                    otherCount += int(record[i])
        dp = refCount + altCount
        if dp > 0:
            vaf = altCount / dp
        else:
            vaf = 0
        res = {'refCount': refCount, 'altCount': altCount, 'otherCount': otherCount,
                    'dp': dp, 'vaf': vaf}
        return  res

    # def get_vaf_dp(self, variantAllele):
    #     res = dict()
    #     upperHeader = [x.upper() for x in self.header]
    #     for key, record in self.data.iteritems():
    #         refCount, altCount, otherCount = 0, 0, 0
    #         if key in variantAllele:
    #             ref, alt = variantAllele[key]
    #             for i in self.ntIdx:
    #                 if upperHeader[i] == ref:
    #                     refCount += int(record[i])
    #                 elif upperHeader[i] == alt:
    #                     altCount += int(record[i])
    #                 else:
    #                     otherCount += int(record[i])
    #             dp = refCount + altCount
    #             if dp > 0:
    #                 vaf = altCount / dp
    #             else:
    #                 vaf = 0
    #             res[key + '_%s_%s'%(ref, alt)] = {'refCount': refCount, 'altCount': altCount, 'otherCount': otherCount,
    #                         'dp': dp, 'vaf': vaf}
    #     return  res

    def output_vaf_dp(self, variantFn, outputFn):
        variantAllele = get_allele(variantFn)
        upperHeader = [x.upper() for x in self.header]
        with open(outputFn, 'w') as fout:
            fout.write('\t'.join(['chr', 'pos', 'ref', 'alt', 'dp', 'vaf', 'incompatible']) + '\n')
            for record in self.data:
                pos = [record[i] for i in self.keyIdx]
                key = '_'.join(pos)
                refCount, altCount, otherCount = 0, 0, 0
                if key in variantAllele:
                    ref, alt = variantAllele[key]
                    for i in self.ntIdx:
                        if upperHeader[i] == ref:
                            refCount += int(record[i])
                        elif upperHeader[i] == alt:
                            altCount += int(record[i])
                        else:
                            otherCount += int(record[i])
                else:
                    ref, alt = 'NA', 'NA'
                dp = refCount + altCount
                if dp == 0:
                    fout.write('\t'.join([pos[0], pos[1], ref, alt, '%s'%dp, 'NA', '%s'%otherCount]) + '\n')
                else:
                    fout.write('\t'.join([pos[0], pos[1], ref, alt, '%s'%dp, '%f'%(altCount / dp), '%s'%otherCount]) + '\n')



if __name__ == '__main__':
    import  argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-pileup', type=str, required=True, help='pileup table, from pileupCount.py')
    parser.add_argument('-indel', type=str, default=None, help='indel count table, from abmIndelCount.py')
    parser.add_argument('-variant', type=str, required=True, help='variant file, vcf or maf')
    parser.add_argument('-output', type=str, required=True, help='output file')

    args = vars(parser.parse_args())

    if args['indel'] is None:
        pileup = PileupTable(args['pileup'])
        pileup.output_vaf_dp(args['variant'], args['output'])
    else:
        maf = MAF(args['variant'])
        pileup = PileupTable(args['pileup'])
        indel = IndelCountTable(args['indel'])
        # vafTb = pileup.get_vaf_dp(maf.get_ref_alt())
        # vafTb.update(indel.get_vaf_dp())
        maf.get_vaf_table(pileup, indel, args['output'])


