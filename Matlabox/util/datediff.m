function delta = datediff(date1, date2, dateformat)

    if ischar(date1) && ischar(date2)
        delta = datenum(date2, dateformat) - datenum(date1, dateformat);
    else
        if ischar(date1) 
            date1 = {date1};
        end
        if ischar(date2)
            date2 = {date2};
        end
        if ~iscolumn(date1)
            date1 = data1';
        end
        if ~iscolumn(date2)
            date2 = data2';
        end
        
        delta = bsxfun(@minus, datenum(date2, dateformat), datenum(date1, dateformat) );
    end