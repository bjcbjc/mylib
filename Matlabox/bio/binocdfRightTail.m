function pval = binocdfRightTail(x, N, p, varargin)

    upper = false;
    if ~isempty(varargin)
        if strcmp(varargin{1}, 'upper')
            upper = true;
        end
    end

    x = x(:);
    n = length(x);
    if length(N) < n
        N = repmat(N, n, 1);
    end
    if length(p) < n
        p = repmat(p, n, 1);
    end
    
    right = x./N > p;
    pval = NaN(size(x));
    pval(right) = binocdf(x(right), N(right), p(right), 'upper');
    pval(~right) = 1 - binocdf(x(~right)-1, N(~right), p(~right));
    
    
    