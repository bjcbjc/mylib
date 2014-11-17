
ENSEMBL2ENTREZ_HUGO <- function(ens, dictfn, entrez_filter=vector()) {
    dict = read.table(dictfn, header=T, sep='\t', colClasses='character')
    if (length(entrez_filter) > 0) {
        dict = dict[ dict[, 'EntrezGene.ID'] %in% entrez_filter, ]
    }
    idx = match(ens[,1], dict[,'Ensembl.Gene.ID'])
    ens[,2] = dict[idx, 'EntrezGene.ID']
    ens[,3] = dict[idx, 'HGNC.symbol']
    ens[is.na(ens[,2]),2] = ''
    ens[is.na(ens[,3]),3] = ''
    colnames(ens) = c('ENSEMBL', 'ENTREZID', 'SYMBOL')
    return(ens) 
}
