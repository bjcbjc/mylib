function [u, m, n] = uniqueCellRows(cellarrays)
    %because matlab ignores 'rows' in unique() when input is cell array
    %
    
    n = size(cellarrays,2);
    txt = cellarrays(:,1);
    for i = 2:n
        txt = strcat(txt, cellarrays(:,i));
    end
    
    [~, m, n] = unique(txt);
    u = cellarrays(m,:);