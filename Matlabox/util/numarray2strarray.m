function strarray = numarray2strarray(numarray)
    [m, n] = size(numarray);
    if m == 1 && n > 1
        singlecol = true;
    elseif m > 1 && n == 1
        numarray = numarray'; %make it row
        singlecol = true;
    elseif m > 1 && n > 1 %multiple columns
        singlecol = false;
    end
    
    if singlecol
        strarray = textscan( num2str(numarray), '%s');
        strarray = strarray{1};
    else        
        if m>n
            numarray = numarray';            
        end        
        [m1, n1] = size(numarray);            
        strarray = cell(m1,n1);
        for i = 1:m1
            tmp = textscan( num2str(numarray(i,:)), '%s');
            strarray(i,:) = tmp{1};
        end
        if m1 ~= m
            strarray = strarray';
        end
    end