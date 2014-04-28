classdef RNACOUNT < handle
    
    methods (Static)
        function normalizedCount = getDESeqNormalizedCount(cohortCount, sampleCount)
            if nargin < 2, sampleCount = cohortCount; end
            normalizedCount = bsxfun(@rdivide, sampleCount, RNACOUNT.getDESeqSizeFactor(cohortCount, sampleCount));
        end
        function sf = getDESeqSizeFactor(cohortCount, sampleCount)
            if nargin < 2, sampleCount = cohortCount; end
            ref = RNACOUNT.getDESeqReferenceCount(cohortCount);
            ratio = bsxfun(@rdivide, sampleCount, ref);
            ratio(isinf(ratio) ) = NaN;
            sf = nanmedian( ratio, 1);
        end
        function refCount = getDESeqReferenceCount(cohortCount)
            refCount = geomean( cohortCount, 2);
        end
        function tpm = getTPM(count, geneLength, pseudoCount)
            %count: gene x sample
            %geneLength: gene x 1            
            if nargin < 3, pseudoCount = 0.9; end
            vl = bsxfun(@rdivide, bsxfun(@rdivide, count+pseudoCount, sum(count,1)), geneLength);
            tpm = bsxfun(@rdivide, vl, sum(vl, 1) ) * 1e6;
        end
    end
end