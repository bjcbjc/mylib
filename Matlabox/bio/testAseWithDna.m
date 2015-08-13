function res = testAseWithDna(rnaReadCount, dnaReadCount, varargin)
    % rnaReadCount: nvariant x 2; (ref, alt)
    % dnaReadCount: nvariant x 2; (ref, alt)
    % method: {'vafDna', 'vafDiff'}
    %   vafDna: binomial, with dna vaf as p; vafDiff: normal, diff btw vaf of dna and rna
    % direction: 1, -1 or 0
    %
    %
    
    para.expectedAEF = [];
    para.minRnaRead = 6;    
    para.method = {'vafDna', 'vafDiff'};  
    
    para = assignpara(para, varargin{:});
    
    if ischar(para.method)
        para.method = { para.method };
    end
    
    [nSNV, ncol] = size(rnaReadCount);
    if ncol ~= 2
        error('rnaReadCount must have two columns: (ref, alt)');
    end
    if size(dnaReadCount) ~= size(rnaReadCount)
        error('rnaReadCount and dnaReadCount must be of the same size');
    end
    
    rnaTotalReadCount = sum(rnaReadCount, 2);
    dnaTotalReadCount = sum(dnaReadCount, 2);
    
    rnaVaf = rnaReadCount(:,2) ./ rnaTotalReadCount;
    dnaVaf = dnaReadCount(:,2) ./ dnaTotalReadCount;
    
    if ismember('vafDna', para.method)                        
        res.vafDnaPval = NaN( nSNV, 1);
        res.vafDnaPvalRefOverExp = NaN( nSNV, 1);
        res.vafDnaPvalAltOverExp = NaN( nSNV, 1);
        valid = rnaTotalReadCount >= para.minRnaRead;
        res.vafDnaPvalAltOverExp(valid) = binocdfRightTail(rnaReadCount(valid,2), rnaTotalReadCount(valid), dnaVaf(valid));
        res.vafDnaPvalRefOverExp(valid) = 1 - res.vafDnaPvalAltOverExp(valid);
%         res.vafDnaPvalAltOverExp(valid) = 1 - binocdf(rnaReadCount(valid,2)-1, rnaTotalReadCount(valid), dnaVaf(valid));
%         res.vafDnaPvalRefOverExp(valid) = 1 - binocdf(rnaReadCount(valid,1)-1, rnaTotalReadCount(valid), 1-dnaVaf(valid));
        res.vafDnaPval = min(res.vafDnaPvalAltOverExp, res.vafDnaPvalRefOverExp);
    end
    if ismember('vafDiff', para.method)
        res.vafDiff = rnaVaf - dnaVaf;
        vafRna2 = (rnaReadCount(:,2) + 0.5) ./ (rnaTotalReadCount + 1);
        vafDna2 = (dnaReadCount(:,2) + 0.5) ./ (dnaTotalReadCount + 1);
        v = vafRna2 .* (1-vafRna2) ./ (rnaTotalReadCount + 1) + ...
            vafDna2 .* (1-vafDna2) ./ (dnaTotalReadCount + 1);                
        res.vafDiffPval = 2 * normcdf(-abs(res.vafDiff), 0, sqrt(v));
        res.vafDiffPvalAltOverExp = normcdf(res.vafDiff, 0, sqrt(v), 'upper');
        res.vafDiffPvalRefOverExp = normcdf(res.vafDiff, 0, sqrt(v));
%         res.vafDiffPval = normcdf(res.vafDiff, 0, sqrt(v), 'upper');
    end
    
    
    