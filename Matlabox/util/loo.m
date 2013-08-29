function [yhat bhat] = loo(y, X)
    
    vi = find(~isnan(y));
    n = length(vi);
    
    yhat = NaN(size(y));
    bhat = NaN(size(X,2), n);
    
    for i = 1:n
        j = vi(setdiff(1:n,i));
        bhat(:,i) = regress(y(j), X(j,:));
        yhat(vi(i)) = X(vi(i),:) * bhat(:,i);
    end
    