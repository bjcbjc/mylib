
import gtf
import argparse

parser = argparse.ArgumentParser();
parser.add_argument('-gtf', type=str, required=True, help='GTF file')
parser.add_argument('-direction', type=str, default='3', choices=['3','5'], help='3 or 5')
parser.add_argument('-skipShorter', type=bool, default=True, help='Skip features shorter than desired length')
parser.add_argument('-length', type=int, required=True, help='desired length')
parser.add_argument('-output', type=str, required=True, help='output file name')

args = vars( parser.parse_args())
    
gtfFile = gtf.GTF(args['gtf'])

with open(args['output'], 'w') as fout:
    fout.write('GeneID\tChr\tStart\tEnd\tStrand\n')
    for gene in gtfFile.getGenes():
        newRegion = gene.getEndRegion(args['length'], args['direction'], args['skipShorter'])
        for x, y in newRegion:
            fout.write('%s\t%s\t%d\t%d\t%s\n'%(gene.id, gene.chrm, x, y, gene.strand))
            
    