import collections
__author__ = 'bjchen'


class FusionCatcherOutput(object):
    def __init__(self, fn, columnMapper = None):
        """

        :param fn: file name, fusion catcher's output
        :return:
        """
        self.fn = fn
        self.records = []
        with open(fn, 'rU') as f:
            self.header = f.readline().strip().split('\t')
            for line in f:
                self.records.append(line.strip().split('\t'))
        self.idxtb = { k:idx for idx, k in enumerate(self.header)}
        self.columnMapper = self.set_column_mapper(columnMapper)


    @staticmethod
    def set_column_mapper(columnMapper):
        if columnMapper is None:
            columnMapper = {'gene1': 'Gene_1_symbol(5end_fusion_partner)',
                            'gene2': 'Gene_2_symbol(3end_fusion_partner)',
                            'loc5': 'Fusion_point_for_gene_1(5end_fusion_partner)',
                            'loc3': 'Fusion_point_for_gene_2(3end_fusion_partner)',
                            'seq': 'Fusion_sequence',
                            'support': 'Spanning_pairs'}
        return columnMapper

    @staticmethod
    def get_id(counter, gene):
        counter[gene] += 1
        recordId = '%s_fus.%d'%(gene, counter[gene])
        return recordId

    @staticmethod
    def get_base(seq, strand5, strand3):
        comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
        seq = seq.split('*')
        base5, base3 = seq[0][-1], seq[1][0]
        if strand5 == '-':
            base5 = comp[base5]
        if strand3 == '-':
            base3 = comp[base3]
        return base5, base3


    def get_data(self, record, key):
        if type(key) is str:
            return record[self.idxtb[key]]
        elif type(key) is list or type(key) is set:
            return [record[self.idxtb[k]] for k in key]

    def convert_to_vcf_record(self, outFn, sampleName):
        # we need, 5 position/strand/base, 3 position/strand/base
        # qual?
        # ID?
        # filter = PASS
        # info: svtype=BND;MATEID=ID?
        ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
        ##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
        ##FORMAT=<ID=SPP,Number=1,Type=Integer,Description="Number of paired-end reads supporting the fusion">

        idCounter = collections.Counter()

        keys = [self.columnMapper[x] for x in ['gene1', 'gene2', 'loc5', 'loc3', 'seq', 'support']]
        with open(outFn, 'w') as fout:
            fout.write('##fileformat=VCFv4.1\n')
            fout.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
            fout.write('##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">\n')
            fout.write('##FORMAT=<ID=SPP,Number=1,Type=Integer,Description="Number of paired-end reads supporting the fusion from fusion catcher">\n')
            fout.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sampleName]) + '\n')
            for record in self.records:
                gene1, gene2, loc5, loc3, seq, support = self.get_data(record, keys)
                chrm5, pos5, strand5 = loc5.split(':')
                chrm3, pos3, strand3 = loc3.split(':')
                base5, base3 = self.get_base(seq, strand5, strand3)

                if strand5 == strand3: # 5[3[ and ]5]3
                    alt5 = '{base5}[{chrm3}:{pos3}['
                    alt3 = ']{chrm5}:{pos5}]{base3}'
                elif strand5 == '+':  # + / -; 5]3] and 3]5]
                    alt5 = '{base5}]{chrm3}:{pos3}]'
                    alt3 = '{base3}]{chrm5}:{pos5}]'
                else: # - / +; [3[5 and [5[3
                    alt5 = '[{chrm3}:{pos3}[{base5}'
                    alt3 = '[{chrm5}:{pos5}[{base3}'

                alt5 = alt5.format(base5=base5, chrm3=chrm3, pos3=pos3)
                alt3 = alt3.format(base3=base3, chrm5=chrm5, pos5=pos5)
                id5 = self.get_id(idCounter, gene1)
                id3 = self.get_id(idCounter, gene2)

                fout.write('\t'.join([chrm5, pos5, id5, base5, alt5, '.', 'PASS', 'SVTYPE=BND;MATEID=%s'%id3, 'SPP', support]) + '\n')
                fout.write('\t'.join([chrm3, pos3, id3, base3, alt3, '.', 'PASS', 'SVTYPE=BND;MATEID=%s'%id5, 'SPP', support]) + '\n')
        return



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-input', type=str, required=True, help='fusion result')
    parser.add_argument('-output', type=str, required=True, help='output VCF')
    parser.add_argument('-sampleName', type=str, required=True, help='sample name')

    args = vars(parser.parse_args())

    fc = FusionCatcherOutput(args['input'])
    fc.convert_to_vcf_record(args['output'], args['sampleName'])