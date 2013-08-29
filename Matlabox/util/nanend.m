function endcol = nanend(mtx)
    %find the index of the last non-NaN column
    %will be useful when sequentially filling a matrix initialized as NaN
    
    [n1 n2] = size(mtx);
    
    if n2 > 1
        endcol = find(isnan(mtx(1,:)), 1, 'first') - 1;
    else
        endcol = find(isnan(mtx), 1, 'first') - 1;
    end

