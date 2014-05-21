function [sortedData, sortIdx] = sortBy(data, label)

    if min(size(data)) > 1
        error('data should be a vector');
    end
    if ~all(ismember(label, data)) || ~all(ismember(data, label))
        error('data and label do not overlap');
    end
    [u, ~, uidx] = unique(data);
    [~, si] = ismember(label, u);
    
    sortedData = data;
    sortIdx = NaN(size(data));
    curIdx = 0;
    for i = 1:length(label)
        idx = uidx == si(i);
        sortIdx(curIdx+1:curIdx+sum(idx)) = find(idx);
        sortedData(curIdx+1:curIdx+sum(idx)) = label(i);
        curIdx = curIdx + sum(idx);
    end
    