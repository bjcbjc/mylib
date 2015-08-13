modzscore = function(observedData, cohortData, MARGIN) {
    # R's mad already applies 1.486 factor
    medAbsDev = apply(cohortData, MARGIN, FUN= function(x) mad(x, na.rm= TRUE))
    meanAbsDev = apply(cohortData, MARGIN, FUN= function(x) mad(x, center= mean(x), constant= 1.253314, na.rm= TRUE))
    idx = medAbsDev == 0
    medAbsDev[idx] = meanAbsDev[idx]

    m = apply(cohortData, MARGIN, FUN= function(x) median(x, na.rm= TRUE))
    tmp = sweep(observedData, MARGIN, m, FUN= '-')
    modz = sweep(tmp, MARGIN, medAbsDev, FUN= '/')

    if (any(medAbsDev == 0) ) {
        idx = sweep(tmp > 0, MARGIN, medAbsDev == 0, FUN= '&')
        modz[idx] = Inf
        idx = sweep(tmp < 0, MARGIN, medAbsDev == 0, FUN= '&')
        modz[idx] = -Inf
    }
    return (modz)
}
