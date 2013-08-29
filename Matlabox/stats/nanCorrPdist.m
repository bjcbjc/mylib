function d = nanCorrPdist(vec, mtx)

    d = 1 - corr(mtx', vec', 'rows', 'pairwise');
    d(isnan(d)) = -2;