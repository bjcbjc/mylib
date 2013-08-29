function pval = lookupPvalTable(stats, pvalTableBin, pvalTable)
    %
    %stats: vectors of statistics used to find p-values in pvalTable
    %pvalTableBin: vector corresponding to the statistics in pvalTable
    %pvalTable: vector of p-values
    %

    if size(stats, 1) == 1
        stats = stats';
    end
    
    if size(pvalTableBin, 1) == 1
        pvalTableBin = pvalTableBin';
    end
        
    n = length(stats);
    
    if ~issorted(pvalTableBin)
        [pvalTableBin, si] = sort(pvalTableBin);
        pvalTable = pvalTable(si);
    end
    
    pval = NaN(n, 1);
    [~, si] = sort( [ stats; pvalTableBin ] );
    binidx = find(si > n);    
    for i = 1:length(binidx)-1
        idx = binidx(i)+1 : binidx(i+1)-1;
        if isempty(idx), continue; end
        pval( si( idx ) ) = pvalTable(i);
    end
    pval( si( 1 : binidx(1)-1 ) ) = pvalTable(1);
    pval( si( binidx(end)+1 : end ) ) = pvalTable(end);
end