function sd = trimstd(X, percent, flag, dim)
    
    if nargin < 4, dim = 1; end
    if nargin < 3, flag = 0; end
    
    if any(isnan(X(:)))        
        error('NaN in the data');        
    end
    
    X = sort(X, dim);
    n = size(X, dim);
    ex = round( n * percent/100/2 );
    
    if n-2*ex == 0
        error(n-ex-ex==0, 'no data left');
    elseif n-2*ex < 3
        warning('<3 data point');
    end
    
    if dim == 1
        sd = std(X( ex+1:n-ex, :), flag, dim);
    else
        sd = std(X( :, ex+1:n-ex), flag, dim);
    end
        