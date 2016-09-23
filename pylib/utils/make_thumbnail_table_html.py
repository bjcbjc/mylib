
import argparse
import glob
import itertools

parser = argparse.ArgumentParser(argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-filePattern', type=str, required=True, help='paths to image files')
parser.add_argument('-nCol', type=int, default=8, help='number of columns for the "image table"')
parser.add_argument('-output', type=str, required=True, help='output file')
parser.add_argument('-height', type=int, default=145, help='height of image')
parser.add_argument('-width', type=int, default=145, help='width of image')

args = vars(parser.parse_args())

img = '<td><A href={src}><img height={height} width={width} src={src}></A></td>'
fns = glob.glob(args['filePattern'])

n = len(fns)
with open(args['output'], 'w') as f:
    f.write('<html><table>\n')
    for i in itertools.islice(xrange(n), 0, n+1, args['nCol']):
        f.write('<tr>\n')
        for j in itertools.islice(fns, i, i+args['nCol']):
            f.write(img.format(src=j, **args) + '\n')
        f.write('</tr>\n')
    f.write('</table></html>\n')
