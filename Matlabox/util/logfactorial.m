function l = logfactorial(k, base)
    %k: scalar or vector
    %calculate log( factorial )
    %
    
    
    if nargin < 2
        d = 1;
    else
        d = log(base);
    end
    
    if length(k) > 1
        K = max(k);
        c = cumsum(log(1:K));
        [~,i] = ismember(k, 1:K);
        l = c(i);
    else
        l = sum(log(1:k));
    end
    
    l = l./d;
    