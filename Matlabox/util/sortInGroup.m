function [sortdata, sortidx] = sortInGroup(data, group, direction)
    if nargin < 3, direction = 'ascend'; end
    
    [n1, n2] = size(data);
    
    if n1 > 1 && n2 > 1
        error('data must be a vector');
    end
    
    if any( sort([n1, n2]) ~= sort(size(group)) )
        error('data and group must have the same dimension');
    end
    
    [~, ~, ui] = unique(group);
    sortidx = zeros( max([n1, n2]), 1);
    cur = 0;
    for i = 1:length(ui)
        si = find(ui == i);
        [~, ii] = sort(data(si), direction);
        groupsize = length(si);
        sortidx(cur+1 : cur+groupsize) = si(ii);
        cur = cur + groupsize;
    end
    
    sortdata = data(sortidx);
    