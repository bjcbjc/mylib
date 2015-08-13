function res = testAseWithPurityEstimate(altReadCount, totalReadCount, varargin)
    %altReadCount: nvariant x 1
    %totalReadCount: nvariant x 1
    %binop: parameter for binomial
    % direction: %1: over-exp of alt, -1: over-exp of ref, 0: both    
    %   1 or -1: one-tail test for over-expression of the allele
    
    para.expectedAEF = [];
    para.minAltRead = 6;
    para.minRead = 6;
    para.minAEF = 0.02;
    para.minN = 30;
    
    
    para = assignpara(para, varargin{:});
    
    aef = altReadCount ./ totalReadCount;
    nloc = length(aef);
%     valid = altReadCount >= para.minAltRead | aef > para.minAEF;
    valid = altReadCount >= para.minAltRead | min(1-aef, aef) > para.minAEF;
    
    if sum(valid) < para.minN
        res.expectedAEF = NaN;
        res.estimatedAEF = false;
        res.purityPvalAltOverExp = NaN( nloc, 1);
        res.purityPvalRefOverExp = NaN( nloc, 1);
        res.purityPval = NaN( nloc, 1);
        return
    end
    
    if isempty(para.expectedAEF)
        [density, xi] = ksdensity(aef(valid));
        %find mode        
        res.expectedAEF = xi( density == max(density) );
        res.estimatedAEF = true;
        if res.expectedAEF > 0.5;
            res.expectedAEF = 0.5;
            res.estimatedAEF = false;
        end
    else
        res.expectedAEF = para.expectedAEF;
        res.estimatedAEF = false;
    end
    
    valid = totalReadCount >= para.minRead;    
    res.purityPvalAltOverExp = NaN( nloc, 1);
    res.purityPvalRefOverExp = NaN( nloc, 1);
    res.purityPval = NaN( nloc, 1);
    if length(res.expectedAEF) == 1
        res.purityPvalAltOverExp(valid) = binocdfRightTail(altReadCount(valid), totalReadCount(valid), res.expectedAEF);
        res.purityPvalRefOverExp(valid) = 1 - res.purityPvalAltOverExp(valid);
%         res.purityPvalRefOverExp(valid) = binocdfRightTail(totalReadCount(valid) - altReadCount(valid), totalReadCount(valid), 1-res.expectedAEF);
        res.purityPval = min(res.purityPvalAltOverExp, res.purityPvalRefOverExp);
    else
        res.purityPvalAltOverExp(valid) = binocdfRightTail(altReadCount(valid), totalReadCount(valid), res.expectedAEF(valid));
        res.purityPvalRefOverExp(valid) = 1 - res.purityPvalAltOverExp(valid);
%         res.purityPvalAltOverExp(valid) = 1 - binocdf(altReadCount(valid)-1, totalReadCount(valid), res.expectedAEF(valid));
%         res.purityPvalRefOverExp(valid) = 1 - binocdf(totalReadCount(valid) - altReadCount(valid)-1, totalReadCount(valid), 1-res.expectedAEF(valid));
        res.purityPval = min(res.purityPvalAltOverExp, res.purityPvalRefOverExp);
    end
    
    
    
    