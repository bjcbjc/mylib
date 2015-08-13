
rnaUtils.getDESeqNormalizedCount = function(cohortCount, sampleCount= NULL) {
    if (is.null(sampleCount))
        sampleCount = cohortCount
    sizeFactor = rnaUtils.getDESeqSizeFactor(cohortCount, sampleCount)
    sweep(sampleCount, 2, sizeFactor, FUN= '/')
}
rnaUtils.getDESeqSizeFactor = function(cohortCount, sampleCount= NULL) {
    if (is.null(sampleCount))
        sampleCount = cohortCount
    ref = rnaUtils.getDESeqReferenceCount(cohortCount)
    ratio = sweep(sampleCount, 1, ref, FUN= '/')
    ratio[is.infinite(as.matrix(ratio))] = NA
    apply(ratio, 2, function(x) median(x, na.rm= TRUE))
}
rnaUtils.geomean = function(x) {
    if (any(x == 0, na.rm= TRUE))
        res = 0
    else
        res = exp(mean(log(x), na.rm= TRUE))
    return(res)
}
rnaUtils.getDESeqReferenceCount = function(cohortCount) {
    apply(cohortCount, 1, rnaUtils.geomean)
}
                
rnaUtil.getTPM = function(sampleCount, geneLength, pseudoCount = 0) {
    sampleCount = sampleCount + pseudoCount
    data = sweep(sampleCount, 1, geneLength, FUN='/')
    return( sweep(data, 2, colSums(data), FUN='/') * 1e6 )
}
