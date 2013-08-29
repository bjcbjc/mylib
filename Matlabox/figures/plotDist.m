function plotDist(data, npt, varargin)

if nargin < 2, npt = 10; end

[logspec logidx] = ismember('logscale', varargin(1:2:end));
if logspec
    logscale = varargin{ logidx * 2 };
    varargin(logidx*2-1:logidx*2) = [];
else
    logscale = false;
end
    

m = nanmin(data);
M = nanmax(data);
tick = m:(M-m)/npt:M;
[h] = hist(data,tick);

if logscale
    %semilogy(tick(h~=0),h(h~=0)./sum(h), varargin{:});
    plot(tick(h~=0),log10(h(h~=0)./sum(h)), varargin{:});
else
    plot(tick(h~=0),h(h~=0)./sum(h), varargin{:});
end

end
