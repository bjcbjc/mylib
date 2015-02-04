

# First column of samplInfoFile is sample ID, used to match samples in read count (so they should match).
# Usesage: path_rscipt deseq2.template.r _var_=value 
# Example: /R/bin/rscript deseq2.template.r sampleInfoFile=/data/sample.info.txt readCountFile=/data/readocount.txt deseq.design=~condition deseq.varOfInterest=condition deseq.condA=untreated deseq.cond=treated
# Note there cannot be space before and after "="
#
# deseq.varOfInterest, deseq.condA, deseq.condB:
#   If all specified, used for 
#      contrast = c(deseq.varOfInterest, deseq.condB, deseq.condA) 
#   to return return log(B/A) for value A, B in the variable varOfInterest 
#   (eg. varOfInterest=condition, condA=untreated, condB=treated).
#
#   If deseq.condA or deseq.condB are not specified, used for
#      name = deseq.varOfInterest
#   to return the effect of the variable (which can be continuous).
#
#   If multiple comparisons are desired, specified a vector, eg.
#      deseq.condA = c('untreated', 'untreated', 'condA')
#      deseq.condB = c('condA', 'condB', 'condB')
#   for analyzing log(condA/untreated), log(condB/untreated), and log(condB/condA).
#
# filterSample: (important)
#   By default, DESeq2 takes all samples for analysis and return the test results by specifying contrast. However, if the 
#   variable for test has more than two values, contrasts only test for those samples with matching values. Neverthess the less, 
#   all samples are used to estimate dispersion and etc. This result is different from subsetting the data to samples that have 
#   matching values in the contrast before analysis. For example, if condition has three values: untreated, condB, condC. 
#   If the interest is to compare condB to untreated, you can
#      (1) filter out condC samples and use the subset of samples to run DESeq()
#      (2) run DESeq() with all sample and specify contrast = c('condition', 'condB', 'untreated').
#   The results can be different between (1) and (2). Set filterSample to TRUE if (1) is desired.
#
# deseq.test: 'Wald' or 'LRT' (likelihood ratio test)
# deseq.design: GLM formula, eg. ~condition, or ~convariate1+condition
# deseq.reduced: reduced GLM formula, this is for likelihood ratio test only
#
# BJC 07/2014
#


#default variable values

deseq.test='Wald'
deseq.fitType = 'parametric'
deseq.design = '~condition'
deseq.reduced = '~1' #only for LRT; full_formula=deseq.design
deseq.minReplicatesForReplace = 7
deseq.betaPrior = NA #default: true for Wald, false for LRT if not specified
deseq.independentFiltering = TRUE
deseq.alphaForIndependentFiltering = 0.1
deseq.varOfInterest = 'condition'  # variable of interests for comparison
deseq.condA = ''  
deseq.condB = ''

filterSample = FALSE
outputPath = './'
outputTag = ''
sampleInfoFile = ''
readCountFile = ''


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

if ( is.na(deseq.betaPrior) ) {
   deseq.betaPrior = deseq.test == 'Wald'
}

#print out variable values for the record
cat('\n\nVariable values used to the analysis:\n')
for (i in 1:length(defVarName)) {
    cat(defVarName[i], '=', eval(parse(text=defVarName[i])), '\n')
}
cat('\n\n')


###### preapre data
sampleInfo = read.table(sampleInfoFile, sep='\t', header=T, check.names=F) #don't pad X for number-starting column names
readCount = read.table(readCountFile, header=T, sep='\t', row.name=1, check.names=F) #don't pad X for number-starting sample names

#remove NA sample in varOfInterest
if (deseq.test != 'LRT') {
    sampleInfo = sampleInfo[!is.na(sampleInfo[,deseq.varOfInterest]), ]
}

validSample = intersect(colnames(readCount), sampleInfo[,1])
if (length(validSample) == 0) {
   cat('Cannot match samples\n')
   cat('Readcount sample:', colnames(readCount), '\n')
   cat('Sample info:', as.character(sampleInfo[,1]), '\n')
   quit(save='no')
} else {
   cat('\n\nMatched', length(validSample), 'samples in readCountFile and sampleInfoFile\n\n')
}

idx = match(validSample, colnames(readCount))
readCount = readCount[, idx]
idx = match(validSample, sampleInfo[,1])
sampleInfo = sampleInfo[idx, ]

if (any( colnames(readCount) != sampleInfo[,1]) ) {
   cat('Samples not matching\n\n')
   quit(save='no')
}

if (any( readCount %% 1 != 0)) {
   cat('Non integer read count!\n\n')
   cat('Floor the data, but the results might not be valid.\n\n')
   readCount = floor(readCount)
}

if ( length(deseq.condA) != length(deseq.condB) ) {
   cat('deseq.condA and deseq.condB have different numbers of elements\n')
   cat('exit!\n')
   quit(save='no')
}

library(DESeq2)

#### LRT test, all samples
if (deseq.test == 'LRT') { 
    countDataFull = DESeqDataSetFromMatrix(countData=readCount, colData=sampleInfo, design=as.formula(deseq.design))
    countDataFull = DESeq(countDataFull, test=deseq.test, fitType=deseq.fitType, minReplicatesForReplace=deseq.minReplicatesForReplace, betaPrior=deseq.betaPrior, reduced=as.formula(deseq.reduced))
    save(countDataFull, file=paste0(outputPath, outputTag, '.LRT.deseq2.rdata'))

    cat('LRT: ', deseq.design, ' vs ', deseq.reduced, '\n')
    result = results( countDataFull, independentFiltering=deseq.independentFiltering, alpha=deseq.alphaForIndependentFiltering ) 
    comparisonName = paste0('LRT', '_', paste(setdiff(all.vars(as.formula(deseq.design)), all.vars(as.formula(deseq.reduced))), collapse='+'))

    write.table( as.data.frame(result), file=paste0(outputPath, outputTag, '.', comparisonName, '.deseq2.txt'), sep='\t', row.names=T, quote=F)
    png(filename=paste0(outputPath, outputTag, '.', comparisonName, '.MA.png'), type='cairo')
    plotMA(result)
    dev.off()
    png(filename=paste0(outputPath, outputTag, '.', comparisonName, '.dispersion.png'), type='cairo')
    plotDispEsts(countDataFull)
    dev.off()

        
##### subsetting samples in condA and condB
} else if ( filterSample && all(deseq.condA != '') && all(deseq.condB != '')  ) {

    for ( comparisonIdx in 1:length(deseq.condA)) {
     	 condA = deseq.condA[ comparisonIdx ]
     	 condB = deseq.condB[ comparisonIdx ]
         idx = sampleInfo[, deseq.varOfInterest] %in% c(condA, condB)
	 cat('Compare', condB, ' vs ', condA, '\n')
    	 if (sum(idx) < 11) {
             cat('Select samples:', as.character(sampleInfo[idx, 1]), '\n\n')
    	 } else {
    	     cat('Select', sum(idx), 'samples\n\n')
    	 }

	 curSampleInfo = droplevels( sampleInfo[ idx, ] )
	 curReadCount = readCount[, idx] 

	 comparisonName = paste0(deseq.varOfInterest, '_', condB, '_vs_', condA)

	 countDataFull = DESeqDataSetFromMatrix(countData=curReadCount, colData=curSampleInfo, design=as.formula(deseq.design))
    	 countDataFull = DESeq(countDataFull, test=deseq.test, fitType=deseq.fitType, minReplicatesForReplace=deseq.minReplicatesForReplace, betaPrior=deseq.betaPrior, reduced=as.formula(deseq.reduced))

    	 save(countDataFull, file=paste0(outputPath, outputTag, '.', comparisonName, '.deseq2.rdata'))
	 result = results( countDataFull, contrast=c(deseq.varOfInterest, condB, condA), independentFiltering=deseq.independentFiltering, alpha=deseq.alphaForIndependentFiltering ) 	     
	 write.table( as.data.frame(result), file=paste0(outputPath, outputTag, '.', comparisonName, '.deseq2.txt'), sep='\t', row.names=T, quote=F)
         png(filename=paste0(outputPath, outputTag, '.', comparisonName, '.MA.png'), type='cairo')
         plotMA(result)
         dev.off()
         png(filename=paste0(outputPath, outputTag, '.', comparisonName, '.dispersion.png'), type='cairo')
         plotDispEsts(countDataFull)
         dev.off()
    }

##### use all samples
} else {
    countDataFull = DESeqDataSetFromMatrix(countData=readCount, colData=sampleInfo, design=as.formula(deseq.design))
    countDataFull = DESeq(countDataFull, test=deseq.test, fitType=deseq.fitType, minReplicatesForReplace=deseq.minReplicatesForReplace, betaPrior=deseq.betaPrior, reduced=as.formula(deseq.reduced))
    save(countDataFull, file=paste0(outputPath, outputTag, '.deseq2.rdata'))
    for ( comparisonIdx in 1:length(deseq.condA)) {
     	 condA = deseq.condA[ comparisonIdx ]
     	 condB = deseq.condB[ comparisonIdx ]

	 if ( condA != '' && condB != '') {
	     cat('Contrast', condB, ' vs ', condA, '\n')
	     result = results( countDataFull, contrast=c(deseq.varOfInterest, condB, condA), independentFiltering=deseq.independentFiltering, alpha=deseq.alphaForIndependentFiltering ) 
	     comparisonName = paste0(deseq.varOfInterest, '_', condB, '_vs_', condA)
	 } else {
	     result = results( countDataFull, name=deseq.varOfInterest, independentFiltering=deseq.independentFiltering, alpha=deseq.alphaForIndependentFiltering ) 
	     comparisonName = deseq.varOfInterest
	 }

	 write.table( as.data.frame(result), file=paste0(outputPath, outputTag, '.', comparisonName, '.deseq2.txt'), sep='\t', row.names=T, quote=F)
         png(filename=paste0(outputPath, outputTag, '.', comparisonName, '.MA.png'), type='cairo')
         plotMA(result)
         dev.off()
         png(filename=paste0(outputPath, outputTag, '.', comparisonName, '.dispersion.png'), type='cairo')
         plotDispEsts(countDataFull)
         dev.off()

    }
}





