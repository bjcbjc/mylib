function longestRoot = longestRootString(strings)

    strings = strings(:);
    n = cellfun(@length, strings);
    n = min(n);
    if n == 0
        longestRoot = '';
    else
        charArray = cell2mat(cellfun(@(x) x(1:n), strings, 'unif', 0));    
        n = find(all(bsxfun(@eq, charArray, charArray(1,:)), 1), 1, 'last');
        if ~isempty(n)
            longestRoot = strings{1}(1:n);
        else
            longestRoot = '';
        end
    end
    