

enrich.go <- function (testGene, backgroundGene, inputGeneId= 'ENSEMBL', database= 'org.Hs.eg.db', ontologyList = c('BP', 'MF', 'CC'), extractGeneInGO = FALSE, pvalCutoff = 0.01, customMapping = FALSE, customMappingFile = '/nethome/bjchen/DATA/GenomicInfo/ENSEMBL/ensembl_entrez_10302014.txt') {
    library(database, character.only= TRUE)
    library(GOstats)
    eval(parse(text=paste0('databaseObj=',database)))
    eval(parse(text=paste0('goGenes=mappedkeys(',sub('.db', 'GO', database, fixed= TRUE), ')')))

    if ( inputGeneId == 'ENTREZID' ) {
        testGeneVector = testGene
        univGeneVector = backgroundGene
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
            genemap = select(databaseObj, testGene, c('ENTREZID', 'SYMBOL'), inputGeneId)
            univmap = select(databaseObj, backgroundGene, c('ENTREZID','SYMBOL'), inputGeneId)
            testGeneVector = genemap[, 'ENTREZID']
            univGeneVector = univmap[, 'ENTREZID']
        }
    }

    #remove genes not in GO or duplicates
    testGeneVector = unique( testGeneVector[ testGeneVector %in% goGenes ] ) 
    univGeneVector = unique( univGeneVector[ univGeneVector %in% goGenes ] ) 

    enrichRes = list()
    for (ontIdx in 1:length(ontologyList)) {
        param = new('GOHyperGParams', geneIds=testGeneVector, universeGeneIds=univGeneVector, annotation=database, ontology=ontologyList[ontIdx], pvalueCutoff=pvalCutoff, conditional=FALSE, testDirection='over')
        hyp = hyperGTest(param)
        sum_hyp = summary(hyp)
        enrichRes[[ontologyList[ontIdx]]] = sum_hyp
        #write.table(sum_hyp, outputFileName, sep='\t', quote=F, row.names=F)
	
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

            #write.table(binaryTable, paste0(outputTag,'.geneByGo_binaryTable_',ontologyList[ontIdx], '.txt'), sep='\t',quote=F,col.names=NA)
            enrichRes[[paste0(ontologyList[ontIdx], '_binaryTable')]] = binrayTable
        }
    }

    return(enrichRes)
}



