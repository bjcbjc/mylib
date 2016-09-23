import sys
__author__ = 'bjchen'


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-input', type=str, default=None)
    parser.add_argument('-output', type=str, default=None)
    parser.add_argument('-tabDelim', type=bool, default=False)

    args = vars(parser.parse_args())

    if args['input'] is None:
        inStream = sys.stdin
    else:
        inStream = open(args['input'], 'rU')

    if args['output'] is None:
        outStream = sys.stdout
    else:
        outStream = open(args['output'], 'w')

    count = 1
    for line in inStream:
        if args['tabDelim']:
            line = line.lstrip().strip('\n').split('\t')
        else:
            line = line.lstrip().strip('\n').split()
        if count == 1:
            outStream.write('||' + '||'.join(line) + '||\n')
        else:
            outStream.write('|' + '|'.join(line) + '|\n')
        count += 1

    inStream.close()
    outStream.close()