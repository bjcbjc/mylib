function [res] = univar_normalization (x, dim)

if nargin < 2, dim = 2; end
if dim == 2
    %res = x ./ (nanstd(x')'*ones(1,size(x,2)));    
    res = bsxfun(@rdivide, x, nanstd(x,0,2));
else
    %res = x ./ repmat(nanstd(x),size(x,1),1);
    res = bsxfun(@rdivide, x, nanstd(x));
    
end

