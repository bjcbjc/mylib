function res = getMarkerGeneWithSNP(markers, SNP, SGD, window, snptype)
    %get genes in markers. defined by Noel.markergenes
    %snptype = {'csnp','nsnp'...} supported by query SNP
    %
    if nargin < 4, window = 30000; end
    if nargin < 5, snptype = {'csnp','nsnp'}; end
        
    M = length(markers);
    res = cell(M,1);

    nsnptype = length(snptype);
    for i = 1:M
        orfs = geneInMarker(markers{i}, window, 'type','Verified','name','orf','SGD',SGD);
        gi = strIndexQuery(SGD.orf, orfs);
        snps = querySNP(orfs, snptype, SNP);

        mstr = sprintf('%s{',markers{i});
        cflag = false;        
        for j = 1:length(orfs)
            if ~isnan(snps(j,1))
                snpStr = sprintf('(%d',snps(j,1));
            else
                snpStr = '(-';
            end
            for snptypei = 2:nsnptype
                if ~isnan(snps(j,snptypei))
                    snpStr = [snpStr sprintf(',%d',snps(j,snptypei))];
                else
                    snpStr = [snpStr ',-'];
                end
            end
            snpStr = [snpStr ')'];
            if ~cflag
                mstr = [mstr sprintf('%s%s', SGD.gene{gi(j)}, snpStr) ];
                cflag = true;
            else
                mstr = [mstr ',' sprintf('%s%s', SGD.gene{gi(j)}, snpStr) ];
            end
        end
        mstr = [mstr '}'];
        res{i} = mstr;
    end
end
