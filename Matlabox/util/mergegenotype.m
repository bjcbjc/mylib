function [representative tb covered] = mergegenotype(loc, genotype, rthres, distthres)
    %merge genotype
    %rthres is the threshold for corr
    %distthres is the threshold for distance
    %loc should be [chrm, loc1, loc2(if available)]
    %genotype: #marker x #sample
    %
    %representative: indexes that are representatitve for the data
    %tb: representatives and those shared corr with the representatives
    %
    %Greedy search until cover the entire set
    %
    %
    %
    
    if nargin < 4, distthres = Inf; end
    if nargin < 3, rthres = 0.8; end
        
    
    allchrm = unique(loc(:,1));
    if size(loc,2) > 2
        loc(:,2) = nanmean(loc(:,2:3), 2);
        loc(:,3) = [];
    end
    
    nmarker = size(genotype,1);
    covered = false(nmarker,1);
    representative = NaN(nmarker, 1);
    tb = cell(nmarker, 1);
    
    repi = 1;
    for ci = 1:length(allchrm)
        vi = find(loc(:,1) == allchrm(ci));
        r = corr(genotype(vi,:)');
        dist = abs(bsxfun(@minus, loc(vi,2), loc(vi,2)'));
        
        r(r < rthres) = 0; %used to find candidate, keep updating
        r(dist > distthres) = 0;
        
        rcopy = r; %used to find correlated members; a gene can have multiple represetatives
        for i = 1:length(vi)
            rcopy(i,i) = 0;
        end
        
        ncorred = sum( r ~= 0, 2);
        [~, candidate] = max(ncorred); %min(ncorred) >= 1 because of self corr; for convergence
        
        %[~, si] = sort(ncorred, 'descend'); %sorted indexes within chrm
        %sortedvi = vi(si); %sorted indexes to all genotype        
        %sii = 1; %indexes to si and sortedvi
        while ~all(covered(vi)) && nnz(r) > 0%|| sii > length(si)
            if ~covered(vi(candidate))
                covered(vi(candidate)) = true;
                representative(repi) = vi(candidate);
                members = find( rcopy( candidate, :) ~= 0); %all correlated
                covered(vi(members)) = true;
                tb{repi,1} = vi(candidate);
                tb{repi,2} = vi(members);
                
                r( members, :) = 0; r(:, members) = 0;
                repi = repi + 1;
            end
            r( candidate, :) = 0; r(:, candidate) = 0;
            ncorred = sum( r ~= 0, 2);
            [~, candidate] = max(ncorred);
%             if ~covered(sortedvi(sii))
%                 covered(sortedvi(sii)) = true;
%                 representative(repi) = sortedvi(sii);
%                 members = find( r( si(sii), :) ~= 0);
%                 covered( vi(members) ) = true;
%                 tb{repi,1} = sortedvi(sii);
%                 tb{repi,2} = members;
%                 repi = repi + 1;
%             end
%             sii = sii + 1;
        end
    end
    representative(repi:end) = [];
    tb(repi:end, :) = [];
    