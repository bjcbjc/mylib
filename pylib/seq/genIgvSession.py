

from sys import argv

class igv_template:
    resourceStr = '<Resource path="{bam}"/>'
    panelStr = '''<Panel height="135" name="{name}" width="1261">
        <Track altColor="0,0,178" autoScale="true" color="175,175,175" colorScale="ContinuousColorScale;0.0;169.0;255,255,255;175,175,175" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="{bam}_coverage" name="{name} Coverage" showReference="false" snpThreshold="0.2" sortable="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="58.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="{bam}" name="{name}" showSpliceJunctions="false" sortable="true" visible="true">
        </Track>
    </Panel>'''
    content = '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="hg19" version="5">
    <Resources>
        {allResources}
    </Resources>
    {allPanels}
    <Panel height="126" name="FeaturePanel" width="1261">
        <Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" sortable="false" visible="true"/>
        <Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;308.0;255,255,255;0,0,178" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" height="35" id="hg19_genes" name="RefSeq Genes" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="308.0" minimum="0.0" type="LINEAR"/>
        </Track>
    </Panel>
</Session>'''
    def validate(self):
        if len(self.bamlist) != len(self.namelist):
            print 'bamlist and namelist are not of the same length'
            exit()
    def __init__(self, bamlist=[], namelist=[], bamlistfn=''):
        self.bamlist = bamlist
        self.namelist = namelist
        if bamlistfn != '':
            bamlist, namelist = self.readBamList(bamlistfn)
            self.addData(bamlist, namelist)
        self.validate()
    def addData(self, bamlist, namelist):
        self.bamlist.extend(bamlist)
        self.namelist.extend(namelist)
        self.validate()
    def getIgvText(self):
        allResources = []
        allPanels = []
        for bam, name in zip(self.bamlist, self.namelist):
            allResources.append(self.resourceStr.format(bam=bam))
            allPanels.append(self.panelStr.format(name=name, bam=bam))
        return self.content.format(allResources='\n'.join(allResources), allPanels='\n'.join(allPanels))
    def readBamList(self,fn):
        bamlist = []
        namelist = []
        for line in open(fn):
            line = line.split()
            bamlist.append(line[0])
            if len(line) > 1:
                namelist.append(line[1])
            else:
                namelist.append(''.join(line[0].split('/')[-1].split('.bam')[:-1]))
        return bamlist, namelist
            
            
if __name__ == '__main__':
    igv = igv_template(bamlistfn = argv[1])
    igvtext = igv.getIgvText()
    if len(argv) > 2:
        f = open(argv[2], 'w')
        f.write(igvtext)
        f.close()
    else:
        print igvtext
