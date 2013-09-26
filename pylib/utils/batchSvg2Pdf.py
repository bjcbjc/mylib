
from sys import argv
from os import popen, system

patterns = argv[1]
if len(argv) < 3:
    if '/' in patterns:
        outdir = patterns[:patterns.index('/')]
    else:
        outdir = './'
else:
    outdir = argv[2]

if outdir[-1] != '/': outdir = outdir + '/'

f = popen('ls %s'%(patterns))
fns = f.read().strip('\n').split('\n')
f.close()

for fn in fns:
    #print 'python /nethome/bjchen/programs/CairoSVG-1.0/cairosvg.py %s -f pdf -o %s%s'%(fn, outdir, fn.replace('.svg', '.pdf'))
    fn = fn.replace(' ', '\ ')
    system('python /nethome/bjchen/programs/CairoSVG-1.0/cairosvg.py %s -f pdf -o %s%s'%(fn, outdir, fn.replace('.svg', '.pdf'))) 
