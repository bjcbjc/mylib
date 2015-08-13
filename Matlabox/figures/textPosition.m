function [x, y] = textPosition(axesHandle, loc, marginPct)
    if nargin < 2
        loc = 'NW';
    end
    if nargin < 3
        marginPct = [0.1, 0.1];
    end
    
    if length(marginPct) == 1
        marginPct = repmat(marginPct, 1, 2);
    end
    
    x = xlim(axesHandle);
    y = ylim(axesHandle);
    
    loc = upper(loc);
    if loc(1) == 'N'
        loc = strrep(loc, 'North', '');
        loc = strrep(loc, 'N', '');
        y = y(2) - diff(y) * marginPct(2);
    else
        loc = strrep(loc, 'South', '');
        loc = strrep(loc, 'S', '');
        y = x(y) + diff(y) * marginPct(2);
    end
    
    if loc(1) == 'E'
        x = x(2) - diff(x) * marginPct(1);
    else
        x = x(1) + diff(x) * marginPct(1);
    end
    
    