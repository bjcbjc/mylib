function s = cellarray2str(cellarray, delimiter)
    %delimiter, def = ' '
    %
    
    if nargin < 2
        delimiter = ' '; 
    else
        delimiter = sprintf(delimiter);
    end
    
    s = cellarray{1};
    for i = 2:length(cellarray)
        if ischar(cellarray{i})
            s = [s sprintf('%s%s',delimiter,cellarray{i})];
        else
            s = [s sprintf('%s%g',delimiter,cellarray{i})];
        end
    end
end