function [lc, uc, cdf] = sigCorrBound(data, nperm, lcut, twotail)

if nargin < 3, lcut = 0.005; end
if nargin < 4, twotail = true; end

[n, m] = size(data);
randc = zeros(m*m, nperm);

for i = 1:nperm
    randc(:,i) = reshape(corr(data, data(randperm(n),:), 'rows', 'pairwise'), m*m, 1);
end

if twotail    
    randc = abs(randc);    
    cbin = 0:0.1:1;
else
    cbin = -1:0.1:1;
end

pdf = histc(randc(:), cbin);
pdf = pdf ./ sum(pdf);
cdf = cumsum(pdf(end:-1:1));
cdf = cdf(end:-1:1);

lidx = find(cdf <= lcut, 1, 'first' );
lc = cbin( lidx );
uidx = lidx + find(cdf(lidx+1:end) == 0, 1, 'first');
if isempty(uidx)
    uc = 1;
else
    uc = cbin( uidx );
end


