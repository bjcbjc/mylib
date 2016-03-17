import re
__author__ = 'bjchen'


class PicardMetrics(object):
    def __init__(self, fn):
        self.fn = fn
        self.input = list()
        self.data = list()
        for line in open(self.fn):
            line = line.strip()
            if line.startswith('#'):
                if 'INPUT=' in line:
                    self.input.extend(re.findall('INPUT=(\S+)', line))
            elif len(line) > 0:
                self.data.append(line)

        self.sample = set()
        if len(self.input) > 0:
            for fileName in self.input:
                self.sample.update(re.findall('Sample_(\S+)/', fileName))
        self.sample = list(self.sample)
        if len(self.sample) > 1:
            print '%s has more than one sample: %s' % (self.fn, self.sample)
        return

class PicardAlignMetrics(PicardMetrics):
    def __init__(self, fn):
        super(PicardAlignMetrics, self).__init__(fn)
        header = self.data[0].split()[1:]
        self.header = list()
        table = list()
        for i in xrange(1, len(self.data)):
            line = self.data.split()
            cat = line[0]
            self.header.extend([x + '_' + cat for x in header])
            table.extend(line[1:])
        self.data = table



