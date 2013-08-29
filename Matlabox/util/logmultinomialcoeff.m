function l = logmultinomialcoeff(mk, base)
    %mk: vector of mi
    %calculate log( sum(mk)! / mk(1)! / mk(2)! / mk(3)! /.../ mk(end)! )
    %
    
    
    if nargin < 2
        d = 1;
    else
        d = log(base);
    end
    
    mk = sort(mk); 
    N = sum(mk);
    k = length(mk);
    
    l = sum( log( mk(end)+1:N ) ) - (k-1) * sum( log( 1:mk(1) ) );

    for i = 2:k-1
        l = l - (k-i) * sum( log( mk(i-1)+1 : mk(i) ) );
    end
    l = l/d;
    