
from sys import argv

fn = argv[1]
outfn = argv[2]

t = open(fn).read().split('^')[1:]
t = [ line for line in t if line.split(' =')[0] == 'SAMPLE' ]

fields = [ [a.split(' =')[0] for a in line.split('!')] for line in t]

uniqueFields = set(fields[0])
for fd in fields: uniqueFields.union(fd)
uniqueFields = list(uniqueFields)

#parse data
table = {}
for fd in uniqueFields:
    table[fd] = []

for record in t:
    data = [ line.strip('\n').split(' = ') for line in record.split('!') ]
    for fd in uniqueFields:
        fdData = [ line[1] for line in data if line[0] == fd ]
        if len(fdData) == 0:
            table[fd].append('NA')
        else:
            table[fd].append( '#'.join( fdData ) )

outLines = []
# not-sample sepcific info
universalFds = [ k for k in table if len(set(table[k])) == 1]
uniqueFields = list(set(uniqueFields).difference(universalFds))
for fd in universalFds:
    outLines.append( '#%s=%s'%(fd, table[fd][0]) )

for fd in uniqueFields:
    outLines.append( '%s\t%s'%(fd, '\t'.join( table[fd])))

fout = open(outfn, 'w')
for line in outLines:
    fout.write('%s\n'%line)
fout.close()
    
    
