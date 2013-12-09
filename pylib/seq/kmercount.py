
from string import maketrans
from collections import Counter
from sys import argv
from itertools import product

k = 9
tb = maketrans('ACGT', '0123')
counter = Counter()
inputfile = argv[1] #assume fastq

test = {}
with open(inputfile) as f:
    for lineno, line in enumerate(f):
        if lineno%4 == 1:
            line = line.strip()
            for i in xrange(0, len(line)-k+1):
                if line[i:i+k] not in test:
                    test[line[i:i+k]] = 1
                else:
                    test[line[i:i+k]] += 1
            line = line.translate(tb)
            for i in xrange(0, len(line)-k+1):
                try:
                    idx = int(line[i:i+k], 4)
                except:
                    continue
                counter[ idx ] += 1


for kmer in product('ACGT', repeat=k):
    kmerstr = ''.join(kmer)
    idx = int( kmerstr.translate(tb), 4)
    print '%s\t%s'%(kmerstr, counter[idx])
