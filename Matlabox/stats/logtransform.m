function data = logtransform(data, logfun, factorForZero)
    if nargin < 2
        logfun = @log2;
    end
    if nargin < 3
        factorForZero = 0.1;
    end

    valid = data ~= 0;
    m = min(data( valid ));
    
    data(valid) = logfun( data(valid) );
    s = std(data(valid));
    data(~valid) = normrnd(logfun( factorForZero * m ), s + logfun(factorForZero), nnz(~valid), 1);
    