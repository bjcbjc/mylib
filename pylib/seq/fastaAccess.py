

def readFAI(faifile):
    data = [ line.strip().split() for line in open(faifile) ]
    fai = {}
    for line in data:
        fai[ line[0] ] = map(int, line[1:]) #chrm length, offset, nt/line, byte_len/line
    return fai

def getChrmByteRange(fai, chrm):
    if chrm not in fai:
        if 'chr' in chrm: 
            chrm = chrm.replace('chr', '').replace('M', 'MT')
        else:
            chrm = 'chr' + chrm
            chrm = chrm.replace('chrMT', 'chrM')
    if chrm not in fai:
        print 'chrm not in index'
        return None
    else:
        chrmLen, offset, ntLen, ntBlen = fai[chrm]
        nWholeLine = chrmLen/ntLen
        if chrmLen % ntLen == 0:
            remainder = 0
        else:
            remainder = chrmLen - nWholeLine * ntLen + (ntBlen - ntLen)
        total = nWholeLine * ntBlen + remainder
        return offset, total

def getChrm(fai, fastaFile, chrm):
    if type(fastaFile) == type('str'):
        fastaFile = open(fastaFile)
    offset, total = getChrmByteRange(fai, chrm)
    fastaFile.seek(offset)
    seq = fastaFile.read(total).replace('\n', '')
    if len(seq) != fai[chrm][0]:
        seq = seq.replace('\r', '')
    if len(seq) != fai[chrm][0]:
        print 'cannot read the chrm %s correctly, %d nt read, but %d indexed'%(chrm, len(seq), fai[chrm][0])
        exit()
    return seq
