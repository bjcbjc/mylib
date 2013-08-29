function idx = triIdx( mtxsize, low)

    if nargin < 2
        low = true;
    end
    
    if low
        idx = tril( ones(mtxsize), -1) == 1;
    else
        idx = triu( ones(mtxsize), 1) == 1;
    end