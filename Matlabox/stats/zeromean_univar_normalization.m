function [x] = zeromean_univar_normalization (x, dim, largescale)

if nargin < 2, dim = 2; end
if nargin < 3, largescale = false; end

%res = univar_normalization ( zeromean_normalization (x, dim), dim );
if any(isnan(x(:)))
    x = bsxfun(@minus, x, nanmean(x,dim));
    if largescale
        x = bsxfun(@rdivide, x, sqrt(nanfastvar(x, 0, dim)));
    else
        x = bsxfun(@rdivide, x, nanstd(x, 0, dim));
    end
else
    warning('off', 'MATLAB:divideByZero');
    x = bsxfun(@minus, x, mean(x,dim));
    if largescale
        x = bsxfun(@rdivide, x, sqrt(fastvar(x, 0, dim)));
    else
        x = bsxfun(@rdivide, x, std(x, 0, dim));
    end
    warning('on', 'MATLAB:divideByZero');
end
x(isinf(x)) = NaN;