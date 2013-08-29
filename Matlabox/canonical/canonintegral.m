function [K, H, G, SCOPE] = canonintegral(k, h, g, scope, intscope)
    %canonical intergral
    %
    %k: matrix, d x d, d is the number of variables/scope for this
    %   canonical form
    %h: matrix, d x n, n is the number of observations (or 1 if no
    %observation)
    %g: vector, 1 x n
    %scope: vector of indexes (representing variables), 1 x d; 
    %intscope: vector of indexes represents the variables over which the
    %   intergal is calculated; this vector must be a subset of scope;     
    %
    %k and h should be ordered according to the order of scope; 
    %
    %If g is empty, its calculation is ignored
    %
    %Return:
    %   K: matrix, s x s, s  = d - p, p is the number of variables over
    %       which the integral is calcualted
    %   H: matrix, s x n
    %   G: vector, 1 x n
    %   SCOPE: 1 x s (indexes for variables) 
    %
    
    [d, n] = size(h);
    
    assert(d == size(k, 1) && d == size(k, 2), ...
        'canonintegral: dimension of k and h must match');
    
    %index intscope
    
    assert(length(scope) == d, ...
        'canonintegral: dimesion of k and scope must match');
    assert(~islogical(intscope), ...
        'canointegral: scope contains indexes but intscope is logical');
    y = ismember(scope, intscope); %1 x d, logical
    
    assert(sum(y) == length(intscope), ...
        'canonintegral: intscope must be a subset of scope');
    
    SCOPE = scope(~y); % 1 x s, preserve the relative ordering in scope
    assert(~isempty(SCOPE), 'canonintegral: empty SCOPE');
    
    p = sum(y);
    
    if ~isempty(g)
        G = zeros(1, n);
        v = g + 0.5 * ( p * log(2*pi) - ...
            log(det(k(y, y))) );  % 1 x 1
        for i = 1:n
            G(i) = h(y, i)' * k(y, y) * h(y, i) ; 
        end
        G = 0.5*G + v; %1 x n
    else
        G = [];
    end
    
    x = ~y;
    kinvk = k(x, y) * ( k(y, y) \ eye(p) );
    K = k(x, x) - kinvk * k(y, x);
    H = h(x, :) - kinvk * h(y, :); %s x n
    
    
    
    
    
    
    
    
    