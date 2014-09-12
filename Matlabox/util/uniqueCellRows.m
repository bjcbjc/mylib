function [u, m, n] = uniqueCellRows(cellarrays, sortOrder)
    %because matlab ignores 'rows' in unique() when input is cell array
    %
    if nargin < 2, sortOrder = 'sorted'; end
    
    n = size(cellarrays,2);
    txt = cellarrays(:,1);
    for i = 2:n
        txt = strcat(txt, cellarrays(:,i));
    end
    
    [~, m, n] = unique(txt, sortOrder);
    u = cellarrays(m,:);