function [K, H, G, SCOPE] = canonprod(k, h, g, scope)
    %canonical product
    %
    %k: cell array, each cell contains the k for a canonical form
    %h: cell array, each cell contains the h for a canonical form
    %g: matrix, d x n, d is the number of canonical forms, n is the number
    %   of observations; n = 1 is no obvervations; g can be a d x 1 cell
    %   array as well, if some g is scaler (no observation); but any cell
    %   with a matrix should have the same dimension
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
        'canonprod: dimension of k, h, scope should be the same');
    if ~isempty(g)        
        assert(ncan == size(g, 1), ...
            'canonprod: dimension of g should be the same as k');        
    end
    
    if ~isempty(g)
        if ~iscell(g)
            G = sum(g, 1); %1 x n
        else
            i = cellfun(@isscalar, g);
            G = sum(cell2mat(g(~i)), 1); % 1 x n
            G = G + sum(cell2mat(g(i)), 1); %the latter should be a scalar
        end
    else
        G = [];
    end
    
    n = max(cellfun(@(x) size(x,2), h));
    %get the scope    
    SCOPE = unique(cell2mat(scope')); %1 x s
    s = length(SCOPE);
    K = zeros(s, s);
    H = zeros(s, n);
    for i = 1:ncan
        [~, idx] = ismember(scope{i}, SCOPE);
        K(idx, idx) = K(idx, idx) + k{i};
        %h can be a mixtures of cells with d x 1 or d x n
        H(idx, :) = bsxfun(@plus, H(idx, :),  h{i}); 
    end
    
    
    
    
    
    
    
    
    