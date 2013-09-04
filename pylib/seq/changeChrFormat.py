
from sys import argv, exit

tb = {'chr1':'1','1':'chr1', 'chr2':'2','2':'chr2', 'chr3':'3','3':'chr3', 'chr4':'4','4':'chr4', 'chr5':'5','5':'chr5', 'chr6':'6','6':'chr6', 'chr7':'7','7':'chr7', 'chr8':'8','8':'chr8', 'chr9':'9','9':'chr9', 'chr10':'10','10':'chr10', 'chr11':'11','11':'chr11', 'chr12':'12','12':'chr12', 'chr13':'13','13':'chr13', 'chr14':'14','14':'chr14', 'chr15':'15','15':'chr15', 'chr16':'16','16':'chr16', 'chr17':'17','17':'chr17', 'chr18':'18','18':'chr18', 'chr19':'19','19':'chr19', 'chr20':'20','20':'chr20', 'chr21':'21','21':'chr21', 'chr22':'22','22':'chr22', 'chr23':'23','23':'chr23','X':'chrX','chrX':'X','Y':'chrY','chrY':'Y','MT':'chrM', 'chrM':'MT'}


inputfile = argv[1]
chrformat = argv[2].lower()

f = open(inputfile, 'r')

if chrformat == 'chr':
    appendChr = True
elif chrformat == 'num':
    appendChr = False
else:
    print 'known chrformat %s'%chrformat
    exit()

if '.fa' in inputfile:
    l = f.readline().strip()
    if '>chr' in l and appendChr:
        print 'already >chr format'
        f.close()
        exit()
                
    while l:
        if l[0] == '>':
            if appendChr: #'>' -> '>chr'
                if l[:3] == '>MT':
                    l = '>chrM' + l[3:]
                else:
                    l = l[0] + 'chr' + l[1:]
            else: #'>chr'  -> '>'
                if l[:5] == '>chrM':
                    l = '>MT' + l[5:]
                else:
                    l = l[0] + l[4:]
        print l
        l = f.readline().strip()
    f.close()
elif '.vcf' in inputfile:
    l = f.readline().strip()
    while l:
        if l[0] == '#':
            print l
            l = f.readline().strip()
        else: break
    if l[:3] == 'chr' and appendChr:
        print 'already in chr format'
        f.close()
        exit()

    while l:
        if appendChr:
            if l[:2] == 'MT':
                l = 'chrM' + l[2:]
            else:
                l = 'chr' + l
        else:
            if l[:3] == 'chr':
                if l[:4] == 'chrM':
                    l = 'MT' + l[4:]
                else:
                    l = l[3:]
        print l
        l = f.readline().strip()
    f.close()
else:
    print 'unknown input format'
