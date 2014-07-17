

# First column of samplInfoFile is sample ID, used to match samples in read count (so they should match).
# Usesage: path_rscipt deseq.template.r _var_=value _var_=\'string_value\'   
# Example: /R/bin/rscript deseq.template.r deseq.method=\'blind\' deseq.maxIteration=200 sampleInfoFile=\'/data/sample.info.txt\'
#
# If deseq.formula0 and deseq.formula1 are specified, GLM comparing deseq.formula1 to deseq.formula0 is used to get the effects of 
# the the addition varibles in formula 1. Variables in formula should correspond to column names in sampleInfoFile.
# 
# If only deseq.formula1 is specified, deseq.condA and deseq.condB can be specified for DE between condA and condB in the
# variable specified in the formula. For example, if deseq.formulta1=\'count ~ treatment\' (treatment is a column in sampleInfoFile),
# specify deseq.condA=\'treated\' and deseq.condB='\untreated\' (assume treated and untreated are the labels in the 'treatment'
# column in sampleInfoFile) to test DE between treated and untreated. 
# 
#
#
#
# BJC 07/2014
#


#default variable values
deseq.method = 'pool-CR'
deseq.sharingMode = 'fit-only'
deseq.fitType = 'parametric'
deseq.maxIteration = 2000
deseq.formula1 = 'count ~ condition'
deseq.formula0 = ''
deseq.condA = ''
deseq.condB = ''
outputPath = './'
outputTag = ''
sampleInfoFile = ''
readCountFile = ''


#overwite default variables if passed by command line
defVarName = ls()
args = commandArgs(trailingOnly=T)
for (i in 1:length(args)) {
    if ( unlist(strsplit(args[i], '='))[1] %in% defVarName ) {
        eval(parse(text=args[i]))
    }
}
if ( substr(outputPath, nchar(outputPath), nchar(outputPath)) != '/') { outputPath = paste0(outputPath, '/') }

#print out variable values for the record
cat('\n\nVariable values used to the analysis:\n')
for (i in 1:length(defVarName)) {
    cat(defVarName[i], '=', eval(parse(text=defVarName[i])), '\n')
}
cat('\n\n')

library(DESeq)

sampleInfo = read.table(sampleInfoFile, sep='\t', header=T, check.names=F) #don't pad X for number-starting column names
readcount = read.table(readCountFile, header=T, sep='\t', row.name=1, check.names=F) #don't pad X for number-starting sample names

validSample = intersect(colnames(readcount), sampleInfo[,1])
if (length(validSample) == 0) {
   cat('Cannot match samples\n')
   cat('Readcount sample:', colnames(readcount), '\n')
   cat('Sample info:', as.character(sampleInfo[,1]), '\n')
   quit(save='no')
} else {
   cat('\n\nMatched', length(validSample), 'samples in readCountFile and sampleInfoFile\n\n')
}

idx = match(validSample, colnames(readcount))
readcount = readcount[, idx]
idx = match(validSample, sampleInfo[,1])
sampleInfo = sampleInfo[idx, ]

if (any( colnames(readcount) != sampleInfo[,1]) ) {
   cat('Samples not matching\n\n')
   quit(save='no')
}

if (any( readcount %% 1 != 0)) {
   cat('Non integer read count!\n\n')
   cat('Floor the data, but the results might not be valid.\n\n')
   readcount = floor(readcount)
}

if (deseq.formula0 == '') {
    # simple DE between condA and condB
    testVarName = gsub(' ', '',  unlist(strsplit(deseq.formula1, '~'))[2])
    validSampIdx = sampleInfo[, testVarName] %in% c(deseq.condA, deseq.condB)
    if (sum(validSampIdx) < 11) {
        cat('Select samples:', as.character(sampleInfo[validSampIdx, 1]), '\n\n')
    } else {
    	cat('Select', sum(validSampIdx), 'samples\n\n')
    }
    countDataFull = newCountDataSet(readcount[, validSampIdx], sampleInfo[validSampIdx, testVarName])
    countDataFull = estimateSizeFactors(countDataFull)
    countDataFull = estimateDispersions(countDataFull, method=deseq.method, sharingMode=deseq.sharingMode, fitType=deseq.fitType)

    result = nbinomTest(countDataFull, deseq.condA, deseq.condB)
    save(countDataFull, result, file=paste0(outputPath, outputTag, '.rdata'))
    write.table( result, file=paste0(outputPath, outputTag, '.deseq.txt'), sep='\t', row.names=F, quote=F)

} else { 
    #full model
    countDataFull = newCountDataSet(readcount, sampleInfo)
    countDataFull = estimateSizeFactors(countDataFull)
    countDataFull = estimateDispersions(countDataFull, method=deseq.method, sharingMode=deseq.sharingMode, modelFormula = deseq.formula1, fitType=deseq.fitType)

    glmfit1 = fitNbinomGLMs( countDataFull, as.formula(deseq.formula1), list(maxit=deseq.maxIteration) )
    glmfit0 = fitNbinomGLMs( countDataFull, as.formula(deseq.formula0), list(maxit=deseq.maxIteration) )
    pvalGLM = nbinomGLMTest( glmfit1, glmfit0 )
    padjGLM = p.adjust( pvalGLM, method="BH" )
    glmfit1$pvalGLM <- pvalGLM
    glmfit1$padjGLM <- padjGLM
    save(countDataFull, glmfit1, glmfit0, file=paste0(outputPath, outputTag, '.rdata'))
    write.table( glmfit1, file=paste0(outputPath, outputTag, '.deseq.txt'), sep='\t', row.names=T, quote=F)
}



