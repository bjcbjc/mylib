function [res] = zeromean_normalization (x, dim)
%default dim = 2; each data entry is a row
%if dim = 1; each data entry is a column
if nargin < 2, dim = 2; end
if dim == 2
    %res = (x - nanmean(x')'*ones(1,size(x,2)));
    res = bsxfun(@minus, x, nanmean(x,2)); %faster when big matrix
else
    %res = (x - repmat(nanmean(x),size(x,1),1));
    %res = (x - ones(size(x,1),1)*nanmean(x)); %faster than repmat
    %mx = nanmean(x);
    %res = x - mx(ones(size(x,1),1),:);
    res = bsxfun(@minus, x,  nanmean(x));
end


