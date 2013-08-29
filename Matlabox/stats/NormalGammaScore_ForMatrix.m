function score = NormalGammaScore_ForMatrix(alpha, lambda, data)
%NormalGammaScore
%   score = NormalGammaScore(alpha,lambda, values)
%
%   Computes the normal gamma score over the set of values in matrix
%   "data", where each row is independent
%   using the given alpha and lambda parameters.

sqrs = nansum(data.*data, 2);
count = sum(~isnan(data),2);
s = nansum(data,2);

if(count == 0)
   score = 0;
   return
end

variance = (sum(sqrs,2) - (s.*s./count))./count;
beta = max(1, lambda / (lambda+1) * (alpha-2));
betaplus = beta + variance.*count/2 + ...
        count .* lambda .* (s./count).^2 ./ (2 * (lambda + count));

alphaplus = alpha + count/2;

score = -count .* log(sqrt(2*pi)) + log(lambda ./ (lambda + count))/2 ...
                + gammaln(alphaplus) - gammaln(alpha) ...
                + alpha * log(beta) - alphaplus.*log(betaplus);

