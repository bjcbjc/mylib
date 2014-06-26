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
        function rpkm = getRPKM(count, geneLength, pseudoCount)
            %count: gene x sample
            %geneLength: gene x 1            
            if nargin < 3, pseudoCount = 0.9; end
            vl = bsxfun(@rdivide, bsxfun(@rdivide, count+pseudoCount, sum(count,1)), geneLength);
            rpkm = vl * 1e9;
        end
        function tmm = getTmmNormalizationFacor(count, varargin)
            %scale normalization; Based on Mark Robinson and Gordon Smyth's
            %R code
            para.Mtrim = 0.3;
            para.Atrim = 0.05;
            para.weighted = true;
            para.Acutoff = -1e+10;
            para.refColIdx = [];            
            para = assignpara(para, varargin{:});
            
            count(all(count==0,2),:) = [];
            N = sum(count, 1);
            tmp = bsxfun(@rdivide, count, N);
            if isempty(para.refColIdx) 
                if size(count,2) == 2
                    para.refColIdx = 1;
                else
                    q = quantile(tmp,0.75,1);
                    [~, para.refColIdx] = min(abs(q-mean(q)));
                end
            end            
            
            M = log2( bsxfun(@rdivide, tmp, tmp(:, para.refColIdx)) );
            A = 0.5 * log2( bsxfun(@times, tmp, tmp(:, para.refColIdx)) );            
            tmp = (1-tmp)./count;
            w = 1 ./ bsxfun(@plus, tmp, tmp(:, para.refColIdx));
            
            tmp = bsxfun(@and, count>0, count(:,para.refColIdx)>0) & A > para.Acutoff;
            M(~tmp) = NaN;
            A(~tmp) = NaN;
            rM = tiedrank(M);
            A = tiedrank(A);
            n = sum(tmp,1);
            nMTrim = floor( n.* para.Mtrim);
            nATrim = floor( n.* para.Atrim);
            
            tmp = bsxfun(@le, rM, nMTrim) | bsxfun(@gt, rM, n-nMTrim) ...
                | bsxfun(@le, A, nATrim) | bsxfun(@gt, A, n-nATrim);
            M(tmp) = NaN;
            w(tmp) = NaN;
            if para.weighted
                tmm = 2.^( nansum(M.*w,1) ./ nansum(w,1));
            else
                tmm = 2.^nanmean(M,1);
            end
            tmm = tmm ./ geomean(tmm);
%             M = log2( (count(:,2)./N(2)) ./ (count(:,1)./N(1)) );
%             A = 0.5*(log2(count(:,2)./N(2)) + log2(count(:,1)./N(1)) );
%             var = (N(2)-count(:,2))./N(2)./count(:,2) + ...
%                 (N(1)-count(:,1))./N(1)./count(:,1);
%             mask = count(:,1)>0 & count(:,2)>0 & A>para.Acutoff;
%             M = M(mask);
%             A = A(mask);
%             var = var(mask);
%             n = sum(mask);
%             idxM = [floor(n*para.Mtrim)+1, NaN];
%             idxM(2) = n + 1 - idxM(1);
%             idxA = [floor(n*para.Atrim)+1, NaN];
%             idxA(2) = n + 1 - idxA(1);
%             
%             rM = tiedrank(M);
%             rA = tiedrank(A);
%             mask = rM >=idxM(1) & rM <=idxM(2) & rA >= idxA(1) & rA <= idxA(2);
%             tmm(1) = 1;
%             tmm(2) = 2.^( sum(M(mask)./var(mask)) ./ sum(1./var(mask)));            
        end
        function cpm = getCPM(count, takeLog, libSize, prior)                        
            % based on edgeR R code by McCarthy and Gordon Smyth
            if nargin < 3, libSize = []; end
            if nargin < 4, prior = 0.25; end
            if isempty(libSize)
                libSize = sum(count, 1);
            elseif strcmp(libSize, 'tmm')
                libSize = sum(count, 1);
                libSize = libSize .* RNACOUNT.getTmmNormalizationFacor(count);
            end
            if takeLog
                prior = libSize ./ mean(libSize)*prior;
                libSize = libSize + 2*prior;            
            end
            libSize = libSize * 1e-6;
            if takeLog
                cpm = log2( bsxfun(@rdivide, bsxfun(@plus, count, prior), libSize));
            else
                cpm = bsxfun(@rdivide, count, libSize);
            end
        end
        
    end
end