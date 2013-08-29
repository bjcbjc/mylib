function mtx = linearscale(mtx, scale, dim, ptile)
%linear scale data to scale; ptile is the percentiles used as max and min
%for scaling
%ptile: eg [0.02 0.98]: 2% and 98% percentile as min and max; avoid extreme
%   outlier

if nargin < 4, ptile = []; end %use max and min; 
if nargin < 3, dim = 1; end
if nargin < 2, scale = [-10 10]; end

if isempty(ptile)
    M = nanmax(mtx, [], dim);
    m = nanmin(mtx, [], dim);
else
    M = prctile(mtx, max(ptile), dim);
    m = prctile(mtx, min(ptile), dim);
end

sM = max(scale);
sm = min(scale);

mtx = sm + bsxfun(@rdivide, bsxfun(@minus, mtx, m), M-m) * (sM-sm);