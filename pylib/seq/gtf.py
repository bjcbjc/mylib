import re
from operator import itemgetter

class Transcript(object):    
    def __init__(self, name, strand, chrm, regionPos):
        self.id = name
        self.strand = strand
        self.chrm = chrm        
        self.regionPos = regionPos
    def getEndRegion(self, desiredLen, direction = '3', skipShorter = False):
        #trim the regions to desired length
        newRegionPos = []
        lenToGo = desiredLen
        rev = False        
        if (direction == '3' and self.strand == '-') or (direction == '5' and self.strand == '+'):
            regs = iter(self.regionPos)
        else:
            regs = reversed(self.regionPos)
            rev = True            
        for x, y in regs:            
            segLen = y - x + 1
            if lenToGo < segLen:
                if rev:
                    newRegionPos.append([y-lenToGo+1, y] )
                else:
                    newRegionPos.append([x, x+lenToGo-1])
                lenToGo = 0
            else:
                newRegionPos.append([x,y])
                lenToGo = lenToGo - segLen
            if lenToGo == 0: 
                break
        if lenToGo > 0 and skipShorter:
            return []        
        if rev:
            newRegionPos = newRegionPos[::-1]            
        return newRegionPos
    
class Gene(Transcript):         
    def __init__(self, name, strand, chrm, regionPos):
        super(Gene, self).__init__(name, strand, chrm, regionPos)
        self.regionPos = self.uniqueRegion()
    def uniqueRegion(self):        
        self.regionPos.sort(key=itemgetter(0))
        uniqReg = [ self.regionPos[0] ]
        last = 0
        for x, y in self.regionPos:
            if uniqReg[last][1] >= x: #overlap, merge
                uniqReg[last][1] = max(y, uniqReg[last][1])
            else:
                uniqReg.append([x, y])
                last += 1
        return uniqReg    
    
                
class GTF:    
    trParser = re.compile('transcript_id "([\-\.\w]+)"')
    geneParser = re.compile('gene_id "([\-\.\w]+)"')
    def __init__(self, fn):
        self.fn = fn
    def getTranscripts(self):
        curname, strand, chrm, regPos = '', '', '', []
        lastReturned = curname
        with open(self.fn) as f:
            for line in f:
                if line[0] == '#': continue
                line = line.split('\t')                
                if line[2] != 'exon': continue                
                try:
                    name = self.trParser.findall(line[-1])[0]
                except:
                    print 'cannot find transcript id'
                    exit()
                if curname != name: #new Transcript
                    if curname != '': 
                        #transcript data
                        tr = Transcript(curname, strand, chrm, regPos)
                    else:
                        tr = None
                    #new transcript
                    curname = name
                    strand = line[6]
                    chrm = line[0]
                    regPos = [ map(int, line[3:5]) ]                    
                    if tr != None:                      
                        lastReturned = tr.id                    
                        yield tr          
                else:
                    regPos.append( map(int, line[3:5]) )                    
            #last
            if curname != lastReturned:
                yield Transcript(curname, strand, chrm, regPos)      

    def getGenes(self):
        curname, strand, chrm, regPos = '', '', '', []
        lastReturned = curname
        with open(self.fn) as f:
            for line in f:
                if line[0] == '#': continue
                line = line.split('\t')                
                if line[2] != 'exon': continue                
                try:
                    name = self.geneParser.findall(line[-1])[0]
                except:
                    print 'cannot find gene id'
                    print line
                    exit()
                if curname != name: #new Transcript
                    if curname != '': 
                        #transcript data
                        gene = Gene(curname, strand, chrm, regPos)
                    else:
                        gene = None
                    #new transcript
                    curname = name
                    strand = line[6]
                    chrm = line[0]
                    regPos = [ map(int, line[3:5]) ]                    
                    if gene != None:                      
                        lastReturned = gene.id                    
                        yield gene          
                else:
                    regPos.append( map(int, line[3:5]) )                    
            #last
            if curname != lastReturned:
                yield Gene(curname, strand, chrm, regPos)      
                                

# gtf = GTF('/nethome/bjchen/Projects/bench/genebodycov/pyscript/test/top1k.gtf')
# count = 0
# for tr in gtf.getTranscripts():
#     print tr.id, tr.strand, tr.chrm, tr.regionPos
#     count = count + 1
# print count
          