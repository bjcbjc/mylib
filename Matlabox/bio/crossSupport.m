function [supportPCT, supportCount, nvar] = crossSupport(GT, varargin)
    para.excludeMissing = true;
    para = assignpara(para, varargin{:});
            
    if para.excludeMissing && nnz(GT < 0) > 0
        fprintf('crossSupport: exclude sites with any missing genotype!\n');        
        idx = any(GT<0, 2);
        fprintf('crossSupport: %d -> %d\n', length(idx), length(idx)-sum(idx));
        GT(idx, :) = [];
    end
    
    [nvar, nsample] = size(GT);
    
    nSupportSample = NaN(nvar, nsample);
    for sampIdx = 1:nsample
        nSupportSample(:, sampIdx) = sum( bsxfun(@eq, GT(:,sampIdx), ...
            GT(:, setdiff(1:nsample, sampIdx)) ), 2);
    end
    supportCount = NaN(nsample-2, nsample);
    for val = 1:nsample-2
        supportCount(val,:) = sum(nSupportSample >= val, 1);
    end

    supportPCT = supportCount(end:-1:1,:)./nvar;
    