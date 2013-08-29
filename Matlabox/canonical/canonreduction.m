function [K, H, G, SCOPE] = canonreduction(k, h, g, scope, e, escope)
    %canonical reduction by evidence
    %
    %k: matrix, d x d, d is the number of variables/scope for this
    %   canonical form
    %h: vector, d x 1
    %g: scalar
    %scope: vector of indexes (representing variables), 1 x d; 
    %e: evidence, p x n, p is the number of variables that have observations
    %   n is the number of observatoins
    %escope: vector of indexes, 1 x p; 
    %
    %k and h should be ordered according to the order of scope; e
    %should be ordered according to the order of escope
    %
    %If g is empty, its calculation is ignored
    %
    %Return:
    %   K: matrix, s x s, s  = d - p, p is the number of observed variables 
    %   H: matrix, s x n
    %   G: matrix, 1 x n
    %   SCOPE: 1 x s (indexes for variables) 
    %
    
    
    d = size(k, 1);
    
    assert(d == size(k, 2) && d == length(h), ...
        'canonreduction: dimension of k and h must match');
    
    %index escope
    
    assert(length(scope) == d, ...
        'canonreduction: dimesion of k and scope must match');
    assert(~islogical(escope), ...
        'canointegral: scope contains indexes but escope is logical');
    
    [~, y] = ismember(escope, scope); %1 x d, index
    
    assert(all(y~=0), ...
        'canonreduction: escope must be a subset of scope');
    
    SCOPE = scope(~ismember(scope, escope)); % 1 x s, preserving the relative ordering in scope
    assert(~isempty(SCOPE), 'canonreduction: empty SCOPE');
    
    if ~isempty(g)
        G = g + h(y)' * e; %1 x n
        n = size(e, 2);
        %for loop is necessary; since n maybe large, it's fewer
        %operations than diag(e'*k(y,y)*e)
        for i = 1:n
            G(i) = G(i) - 0.5 * e(:,i)'*k(y, y)*e(:,i); 
        end
    else
        G = [];
    end
        
    [~, x] = ismember(SCOPE, scope);
    K = k(x, x); % s x s, s = d - p; 
    H = bsxfun(@minus, h(x), k(x, y) * e); % s x n
    
    