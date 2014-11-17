


geneList = ''
background = ''
outputPath = './'
outputTag = ''
database = 'org.Hs.eg.db'
extractGeneInGO = FALSE
inputGeneId = 'ENSEMBL'
pvalCutoff = 0.01
customMapping = FALSE
customMappingFile = '/nethome/bjchen/DATA/GenomicInfo/ENSEMBL/ensembl_entrez_10302014.txt'
ontologyList = c('BP', 'MF', 'CC')


#overwite default variables if passed by command line
defVarName = ls()
args = commandArgs(trailingOnly=T)
for (i in 1:length(args)) {
    varVal = unlist(strsplit(args[i], '='))
    varVal[2] = gsub('"|\'|\\s', '', varVal[2])
    if ( varVal[1] %in% defVarName ) {
        if ( suppressWarnings(!is.na(as.numeric(varVal[2]))) || suppressWarnings(!is.na(as.logical(varVal[2]))) ) {
            eval(parse(text=args[i]))
	} else if ( length(grep('c(', varVal[2], fixed=T)) > 0 ) {
	    eval(parse(text=paste0(varVal[1],'=',gsub(',', '\',\'', gsub(')', '\')', gsub('(', '(\'', varVal[2], fixed=T), fixed=T), fixed=T))))
	} else {
	    eval(parse(text=paste0(varVal[1],'="',varVal[2],'"')))
	}
    }
}
if ( substr(outputPath, nchar(outputPath), nchar(outputPath)) != '/') { outputPath = paste0(outputPath, '/') }

#print out variable values for the record
cat('\n\nVariable values used to the analysis:\n')
for (i in 1:length(defVarName)) {
    cat(defVarName[i], '=', eval(parse(text=defVarName[i])), '\n')
}
cat('\n\n')


testGene = read.table(geneList, header=F, as.is=T, colClasses='character');
backgroundGene = read.table(background, header=F, as.is=T, colClasses='character');

library(database, character.only=T)
library(GOstats)
eval(parse(text=paste0('databaseObj=',database)))
eval(parse(text=paste0('goGenes=mappedkeys(',sub('.db', 'GO', database, fixed=T), ')')))

if ( inputGeneId == 'ENTREZID' ) {
   testGeneVector = testGene[,1]
   univGeneVector = backgroundGene[,1]
} else {
   if (customMapping && inputGeneId == 'ENSEMBL') {
      source('/nethome/bjchen/BJLib/rscripts/convertENSEMBL2ENTREZ.r')
      genemap = ENSEMBL2ENTREZ_HUGO(testGene, customMappingFile, goGenes)
      univmap = ENSEMBL2ENTREZ_HUGO(backgroundGene, customMappingFile, goGenes)
      #remove empty
      genemap = genemap[ genemap[, 'ENTREZID'] != '', ]
      univmap = univmap[ univmap[, 'ENTREZID'] != '', ]
      testGeneVector = genemap[, 'ENTREZID']
      univGeneVector = univmap[, 'ENTREZID']
   } else { 
      genemap = select(databaseObj, testGene[,1], c('ENTREZID', 'SYMBOL'), inputGeneId)
      univmap = select(databaseObj, backgroundGene[,1], c('ENTREZID','SYMBOL'), inputGeneId)
      testGeneVector = genemap[, 'ENTREZID']
      univGeneVector = univmap[, 'ENTREZID']
   }
}

#remove genes not in GO or duplicates
testGeneVector = unique( testGeneVector[ testGeneVector %in% goGenes ] ) 
univGeneVector = unique( univGeneVector[ univGeneVector %in% goGenes ] ) 
    
for (ontIdx in 1:length(ontologyList)) {
    param = new('GOHyperGParams', geneIds=testGeneVector, universeGeneIds=univGeneVector, annotation=database, ontology=ontologyList[ontIdx], pvalueCutoff=pvalCutoff, conditional=FALSE, testDirection='over')
    hyp = hyperGTest(param)
    sum_hyp = summary(hyp)
    outputFileName = paste0(outputTag, '.GOstats_', ontologyList[ontIdx], '.txt')
    write.table(sum_hyp, outputFileName, sep='\t', quote=F, row.names=F)
	
    if (extractGeneInGO) {
        #extract genes that are in each annotation
	gn2go = select(databaseObj, geneIds(hyp), c('GOALL') )
    	gn2go = gn2go[ gn2go[,'GOALL'] %in% sum_hyp[,1], ]

    	rowName = unique(gn2go[, 'ENTREZID'])
    	colName = unique(gn2go[, 'GOALL'])
    	rowIdx = match(gn2go[, 'ENTREZID'], rowName)
    	colIdx = match(gn2go[, 'GOALL'], colName)
    	idx = cbind(rowid=as.vector(t(rowIdx)), colid=as.vector(t(colIdx)))
    	binaryTable = matrix(0, nrow=length(rowName), ncol=length(colName), dimnames=list(rowName, colName))
    	binaryTable[idx] = 1
	#convert GO ID to terms
    	idx = match(colnames(binaryTable), sum_hyp[,1])
    	colnames(binaryTable) = sum_hyp[idx,'Term']
	#add ENSEMBL/SYMBOL
	idx = match(rownames(binaryTable), genemap[, 'ENTREZID'])
	geneInfo = genemap[idx, c('ENSEMBL', 'SYMBOL')]
	rownames(geneInfo) = rownames(binaryTable)
	binaryTable = cbind(geneInfo, binaryTable)

    	write.table(binaryTable, paste0(outputTag,'.geneByGo_binaryTable_',ontologyList[ontIdx], '.txt'), sep='\t',quote=F,col.names=NA)
    }
}

