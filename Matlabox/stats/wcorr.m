function wr = wcorr(X, y, w)

    assert(length(w)==length(y), 'len inconst for y and w');
    assert(length(w)==size(X,1), 'len inconst for X and w');
    
    sw = sum(w);
%     swx = sum(bsxfun(@times, X, w));
%     swy = sum(w.*y);
%     wr = (sum(bsxfun(@times, bsxfun(@times, X, y), w)) - swx.*swy/sw) ./ ...
%         sqrt((sum(bsxfun(@times, X.^2, w))-(swx.^2)/sw) .* ...
%         (sum(y.^2.*w)-(swy.^2)/sw));
    

    %faster
    nmx = bsxfun(@minus, X, sum(bsxfun(@times, X, w))/sw);
    nmy = y - sum(y.*w)/sw;
    wr = ((w.*nmy)' * nmx) ./ sqrt( ...
        sum(bsxfun(@times, nmx.^2, w)) .* sum(w.*(nmy.^2)) );
     