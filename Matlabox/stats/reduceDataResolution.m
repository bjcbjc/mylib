function rdata = reduceDataResolution(data, dim, numpt, func)
    %reduce the data resolution by summarizing data per numpt, using func
    %func can be mean or median (@nanmean, or @nanmedian)
    %
    %
    
    if nargin < 2
        dim = 2;
    end
    if nargin < 3
        numpt = 1000;
    end
    if nargin < 4
        func = @nanmean;
    end
    
    if length(size(data)) > 2
        error('only support 2-D matrix\n');
    end
    
    wasrow = false;
    if size(data,1) == 1
        data = data';
        wasrow = true;
    end
    
    rdata = [];
    n = size(data, dim);
    if n <= numpt
        fprintf('too few data points for reducing per %d points\n', numpt);
        return
    end
    
    range = 1:numpt:n;
    if (n-range(end))/numpt < 0.5
        range(end) = n;
    else
        range = [range n];
    end %don't really need the whole vector...only the last two elements
    
    npt = length(range)-1;
    
    ndim = size(data, setdiff(1:2,dim));
    rdata = NaN(npt, ndim);
    
    [n1, n2] = size(data);
    for di = 1:ndim
        if dim == 2
            idx = di:n1:n1*n2;
        else
            idx = (di-1)*n1+1:di*n1;
        end
        tmp = reshape(data( idx(1:range(end-1)-1) ), numpt, npt-1);
        rdata(1:npt-1, di) = func(tmp, 1);
        rdata(npt, di) = func(data( idx(range(end-1):range(end)) ) );
    end
    
    if wasrow
        rdata = rdata';
    end
    