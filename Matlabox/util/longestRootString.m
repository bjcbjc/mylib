function longestRoot = longestRootString(strings)

    strings = strings(:);
    n = cellfun(@length, strings);
    n = min(n);
    if n == 0
        longestRoot = '';
    else
        charArray = cell2mat(cellfun(@(x) x(1:n), strings, 'unif', 0));    
        cmpidx = bsxfun(@eq, charArray, charArray(1,:));
        n0 = find(any(cmpidx == 0, 1), 1, 'first');
        if isempty(n0)
            n = find(all(cmpidx, 1), 1, 'last');
        else
            n = n0-1;
        end
        if ~isempty(n)
            longestRoot = strings{1}(1:n);
        else
            longestRoot = '';
        end
    end
    