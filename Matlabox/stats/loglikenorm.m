function l = loglikenorm(x, mu, sigma2)
    % support multiple mean's and sigma2's
    % calculate log likelihood
    %
    % input:
    %   x: data: n1 x n2
    %   mu: mean, scalar, n1 x n2, n1 x 1, or 1 x n2
    %   sigma: variance, scalar, n1 x n2, n1 x 1, or 1 x n2
    %
    % return:
    %   l: n1 x n2, point estimate of log likelihood
    %
    
    [n1, n2] = size(x);
    [mn1, mn2] = size(mu);
    [sn1, sn2] = size(sigma2);
    
    muscalar = isscalar(mu);
    sigmascalar = isscalar(sigma2);
    
    muok = n1 == mn1 | n2 == mn2 | muscalar;
    sigmaok = n1 == sn1 | n2 == sn2 | sigmascalar;
    assert(muok && sigmaok, 'dim mismatch');
    
    if muscalar && sigmascalar
        l = -0.5 * ( log(2*pi) + log(sigma2) + ((x - mu).^2)./sigma2 );
    elseif muscalar
        l = -0.5 * ( log(2*pi) + bsxfun(@plus, log(sigma2),  ...
            bsxfun(@rdivide, (x - mu).^2, sigma2 ) ) );
    elseif sigmascalar
        l = -0.5 * ( log(2*pi) + log(sigma2) + (bsxfun(@minus, x, mu).^2)./sigma2 );
    else
        l = -0.5 * ( log(2*pi) + bsxfun(@plus, log(sigma2), ...
            bsxfun(@rdivide, bsxfun(@minus, x, mu).^2, sigma2) ) );
    end