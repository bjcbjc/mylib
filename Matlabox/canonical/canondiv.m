function [K, H, G, SCOPE] = canondiv(k, h, g, scope)
    %canonical division
    %
    %k: cell array, each cell contains the k for a canonical form
    %h: cell array, each cell contains the h for a canonical form
    %g: matrix, d x n, d is the number of canonical forms, n is the number
    %   of observations; n = 1 is no obvervations
    %scope: cell array, each cell specify the scope of variables for a
    %   canonical form; can be vector of indexes representing nodes; 
    %   each cell is a ROW vector
    %
    %Dimension: for each cell, k: d x d; h = d x n; scope = 1 x d 
    %k and h should be ordered according to the order of scope; e
    %should be ordered according to the order of escope; n is the number of
    %instances (observations)
    %
    %If g is empty, its calculation is ignored
    %
    %Return:
    %   K: matrix, s x s, s is the number of union of all scopes
    %   H: matrix, s x n
    %   G: vector, 1 x n
    %   SCOPE: 1 x s (indexes for variables)  
    %
    
    
    %check dimension
    ncan = length(k);
    assert(ncan == length(h) && ncan == length(scope), ...
        'canondiv: dimension of k, h, scope should be the same');
    if ~isempty(g)
        assert(ncan == size(g, 1), ...
            'canondiv: dimension of g should be the same as k');
    end
    
    if ~isempty(g)
        G = g(1, :) - sum(g(2:end), 1);
    else
        G = [];
    end
    
    n = size(g, 2);
    
    SCOPE = unique(cell2mat(scope')); %1 x s
    s = length(SCOPE);
    K = zeros(s, s);
    H = zeros(s, n);
    
    [~, idx] = ismember(scope{1}, SCOPE);
    K(idx, idx) = k{1};
    H(idx) = h{1};
    
    for i = 2:ncan
        [~, idx] = ismember(scope{i}, SCOPE);
        K(idx, idx) = K(idx, idx) - k{i};
        H(idx, :) = H(idx, :) - h{i};
    end
    
    