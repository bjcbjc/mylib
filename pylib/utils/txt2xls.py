import collections
import xlsxwriter

__author__ = 'bjchen'
## combine tables into an excel xlsx file


class xlsxFile(object):
    def __init__(self, xlsxFn):
        self.name = xlsxFn
        # self.workbook = xlsxwriter.Workbook(self.name)
        self.workbook = xlsxwriter.Workbook(self.name, {'strings_to_numbers': True})

    def add_table(self, sheetName, fn):
        worksheet = self.workbook.add_worksheet(sheetName)
        idx = 0
        for line in open(fn):
            line = line.strip().split('\t')
            worksheet.write_row(idx, 0, line)
            idx += 1
        return

    def close(self):
        self.workbook.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage='txt2xls.py -o <output> -<sheetName1> <table1.txt> ...')

    parser.add_argument('-o', type=str, required=True, help='output file name')


    args = parser.parse_known_args()

    if len(args[1]) == 0:
        print 'Usage: txt2xls.py -o <output.xlsx> -<sheetName1> <table1.txt> ...'
        raise ValueError('Need to specify input tables.')
    elif len(args[1]) %2 == 1:
        print 'Usage: txt2xls.py -o <output.xlsx> -<sheetName1> <table1.txt> ...'
        raise ValueError('Need to specify -<sheetName> <table> in pairs.')
    else:
        dataDict = collections.OrderedDict()
        for i in xrange(0, len(args[1]), 2):
            dataDict[args[1][i].lstrip('-')] = args[1][i+1]

    args = vars(args[0])


    if not args['o'].endswith('.xlsx'):
        args['o'] += '.xlsx'

    xf = xlsxFile(args['o'])
    for sheetName, fn in dataDict.iteritems():
        xf.add_table(sheetName, fn)

    xf.close()

